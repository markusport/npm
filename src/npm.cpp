/*! \file npm.cpp
* \brief Implementation of the natal philopatry model
*/

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <utility>
#include <random>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <regex>
#include "population.h"
#include "visitors.h"


namespace npm {


  rndutils::xorshift128 thread_local RndEng = rndutils::make_random_engine<>();
  

  const char* Version = "0.2.1";
  const char* mating_name[Mating::MATING_MAX] = { "random", "residency" };
  const char* oplacement_name[oPlacement::OPLACEMENT_MAX] = { "back", "sort" };
  const char* ovote_name[oVote::OVOTE_MAX] = { "ignore", "account" };
  const char* bvote_name[bVote::BVOTE_MAX] = { "ignore", "kin", "despotic", "egalitarian", "hierarchical" };


  //! \brief The Simulation
  class Simulation
  {
  public:
    //! \brief creates the simulation
    explicit Simulation(Parameter const& param);

    //! Model loop
    template <Mating MODE, oPlacement PLACEMENT>
    void run();

  private:
    void log(size_t T);
    void clog(size_t T);
    std::ostream& stream_R_header(std::ostream& os) const;
    std::ostream& stream_mean_alleles(std::ostream& os);
    std::ostream& stream_mean_xy(std::ostream& os);
    std::ostream& stream_groupsize(std::ostream& os);
    std::ostream& stream_mean_groupsize(std::ostream& os);
    std::ostream& stream_alleles(std::ostream& os);
    std::ostream& stream_mranks(std::ostream& os);
    std::ostream& stream_xynR(std::ostream& os);
    std::ostream& stream_takeover_stats(std::ostream& os);
    Parameter param_;
    Population pop_;
    TakeoverStats takeover_stats_;
    TakeoverStats takeover_stats_log_;
    TakeoverStats takeover_stats_clog_;
    std::chrono::high_resolution_clock::time_point t0_;
    std::ofstream of_;
  };


  Simulation::Simulation(Parameter const& param)
  : param_(param),
    pop_(param_),
    takeover_stats_{ 0, 0, 0 },
    takeover_stats_log_{ 0, 0, 0 },
    takeover_stats_clog_{ 0, 0, 0 }
  {
    // prepare R file
    fs::create_directories(param_.offile.parent_path());
    of_.open(param_.offile, std::fstream::out | std::fstream::trunc);
    if (!of_) 
    { // complain about file creation failure
      throw std::runtime_error((std::string("Can't create output file ") + fs::canonical(param_.offile).string()).c_str());
    }
    stream_R_header(of_);
    if (param_.oany)
    {
      if (param_.ot) std::cout << "'Time' ";
      if (param_.og) std::cout << "'group-size (f)' ";
      if (param_.om) std::cout << "'group-size (m)' ";
      if (param_.off) std::cout << "'floater (f)' ";
      if (param_.omf) std::cout << "'floater (m)'  ";
      if (param_.oa) std::cout << "'Alleles' ";
      if (param_.oxy) std::cout << "'x y' ";
      if (param_.oto) std::cout << "'takeover' ";
      std::cout << std::endl;
    }
    t0_ = std::chrono::high_resolution_clock::now();
  }


  template <Mating MODE, oPlacement PLACEMENT>
  void Simulation::run()
  {
    size_t T = 0;
    for (; T < param_.ticks; ++T)
    {
      for (auto& patch : pop_.patches())
      {
        patch.do_reproduction<MODE>(param_, pop_.male_floater());
        patch.do_dispersal<PLACEMENT>(param_, pop_.female_floater(), pop_.male_floater());
        patch.do_survival(param_);
      }
      pop_.shuffle_floater(param_);
      pop_.do_floater_survival(param_);
      takeover_stats_ += pop_.do_colonization<MODE>(param_);
      pop_.visit_all([](Individual& ind){ ++ind.age; });

      log(T);
      clog(T);
    }
    // Epilogue - append npm.R to result file
    auto cwd = fs::current_path();
    std::ifstream ifs((cwd / "npm.R").c_str());
    of_ << '\n' << ifs.rdbuf();
  }


  void Simulation::log(size_t T)
  {
    if ((param_.log && (T % param_.log == 0)) || (T == param_.ticks - 1))
    { // append to R file
      of_ << "T <- cbind(T, " << T << ")\n";
      if (!param_.aloglast || T == param_.ticks - 1)
      {
        auto p = of_.precision();
        of_ << std::fixed << std::setprecision(param_.precision);
        stream_alleles(of_);
        stream_xynR(of_);
        stream_mranks(of_);
        auto flags = std::resetiosflags(std::ios::fixed);
        of_.precision(p);
      }
      stream_groupsize(of_);
      stream_takeover_stats(of_);
      takeover_stats_log_ = takeover_stats_;
      of_ << "fFloater <- cbind(fFloater, " << pop_.female_floater().size() << ")\n";
      of_ << "mFloater <- cbind(mFloater, " << pop_.male_floater().size() << ")\n";
      of_ << std::endl;
    }
  }


  void Simulation::clog(size_t T)
  {
    if (param_.oany && ((param_.clog && (T % param_.clog == 0)) || (T == param_.ticks - 1)))
    { // print to console
      auto prec = std::cout.precision();
      std::cout.precision(4);
      if (param_.ot) std::cout << T << '\t';
      stream_mean_groupsize(std::cout) << "  ";
      if (param_.oa) stream_mean_alleles(std::cout) << "  ";
      if (param_.oxy) stream_mean_xy(std::cout) << "  ";
      if (param_.oto) 
      { 
        auto t = (takeover_stats_ - takeover_stats_clog_) / param_.clog;
        std::cout << t.attempt << ' ' << t.takeover << ' ' << t.walkin; 
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      if (param_.oprof) std::cout << "\t  " << std::chrono::duration<double>(t1 - t0_).count();
      t0_ = t1;
      std::cout << std::endl;
      std::cout.precision(prec);
      takeover_stats_clog_ = takeover_stats_;
    }
  }


  std::ostream& Simulation::stream_R_header(std::ostream& os) const
  {
    os << "# Natal philopatry model result file\n";
    os << "# Version " << Version << '\n';
    os << "path <- '" << fs::absolute(fs::path(param_.offile).remove_filename()).generic_string() << "'\n";
    os << "file <- '" << param_.offile.filename() << "'\n";
    os << "rep <- " << param_.rep << "\n\n";
    os << "# Parameter set\n";
    os << "m <- " << param_.m << '\n';
    os << "m0 <- " << param_.m0 << '\n';
    os << "F0 <- " << param_.F0 << '\n';
    os << "phi <- " << param_.phi << '\n';
    os << "delta <- " << param_.delta << '\n';
    os << "k <- " << param_.k << '\n';
    os << "Alleles <- c(" << param_.alleles[0];
    for (size_t i = 1; i < Loci::MAX_ALLELE; ++i) os << ", " << param_.alleles[i]; 
    os << ")\n";
    os << "Mask <- c(" << param_.mask[0];
    for (size_t i = 1; i < Loci::MAX_ALLELE; ++i) os << ", " << param_.mask[i]; 
    os << ")\n";
    os << "Sb <- " << param_.Sb << '\n';
    os << "Sm <- " << param_.Sm << '\n';
    os << "Sff <- " << param_.Sff << '\n';
    os << "Smf <- " << param_.Smf << '\n';
    os << "Smax <- " << param_.Smax << '\n';
    os << "sigma <- " << param_.sigma << '\n';
    os << "thetaB <- " << param_.thetaB() << '\n';
    os << "thetaM <- " << param_.thetaM() << '\n';
    os << "eps <- " << param_.eps << '\n';
    os << "t0 <- " << param_.t0 << '\n';
    os << "tau <-" << param_.tau << '\n';
    os << "mu <- " << param_.mu << '\n';
    os << "mudist <- " << (std::is_same<mutation_dist, std::cauchy_distribution<>>::value ? "'cauchy'\n" : "'normal'\n");
    os << "gamma <- " << param_.gamma << '\n';
    os << "mode <- '" << mating_name[(int)param_.mode] << "'\n";
    os << "ovote <- '" << ovote_name[(int)param_.ovote] << "'\n";
    os << "bvote <- '" << bvote_name[(int)param_.bvote] << "'\n";
    os << "oplacement <- '" << oplacement_name[(int)param_.oplacement] << "'\n";
    os << "ticks <- " << param_.ticks << '\n';
    os << "log <- " << param_.log << "\n";
    os << "aloglast <- " << (param_.aloglast ? 1 : 0) << "\n\n";
    os << "T <- list()        # Vector of log-times\n\n";
    os << "# inherited alleles and response of the breeders per log\n";
    os << "# Each element in the following lists is a matrix(..., nrow = number alleles)\n";
    os << "allele0 <- list()  # list of first allele at gene loci A0, A1, A2, B0, B1, B2 per individual\n";
    os << "allele1 <- list()  # list of second allele at gene loci A0, A1, A2, B0, B1, B2 per individual\n";
    os << "xynR <- list()     # list of x(n,R) and y(n,R) per individual\n\n";
    os << "mrank <- list()    # rank of the breeders mother at birth\n";
    os << "gs <- list()       # group sizes\n";
    os << "males <- list()    # resident males\n";
    os << "takeover <- list() # {attempted, successful, walk-in}\n";
    os << "fFloater <- list() # number of female floater\n";
    os << "mFloater <- list() # number of male floater\n";
    os << std::endl;
    return os;
  }


  std::ostream& Simulation::stream_mean_alleles(std::ostream& os)
  {
    mean_allele_visitor mav;
    pop_.visit_all(std::ref(mav));
    auto mean = mav.mean();
    for (size_t i = 0; i < Loci::MAX_ALLELE; ++i)
    {
      os << mean[i] << ' ';
    }
    return os;
  }


  std::ostream& Simulation::stream_mean_xy(std::ostream& os)
  {
    mean_behavior_visitor mxyv;
    pop_.visit_patches(std::ref(mxyv));
    auto mean = mxyv.mean();
    os << mean.x << ' ' << mean.y << ' ';
    return os;
  }


  std::ostream& Simulation::stream_groupsize(std::ostream& os)
  {
    os << "gs[[length(gs)+1]] = c(";
    for(auto const& patch : pop_.patches()) os << patch.size() << ',';
    os.seekp(-1,std::ios_base::cur) << ")\n";
    os << "males[[length(males)+1]] = c(";
    for(auto const& patch : pop_.patches()) os << (patch.male() == nullptr ? 0 : 1) << ',';
    os.seekp(-1,std::ios_base::cur) << ")\n";
    return os;
  }


  std::ostream& Simulation::stream_mean_groupsize(std::ostream& os)
  {
    size_t s = 0;
    size_t m = 0;
    for (auto const& patch : pop_.patches()) { s += patch.size(); m += (nullptr == patch.male()) ? 0 : 1; }
    if (param_.og) os << static_cast<double>(s) / pop_.patches().size() << ' '; 
    if (param_.om) os << static_cast<double>(m) / pop_.patches().size() << ' ';
    if (param_.off) os << pop_.female_floater().size() << ' ';
    if (param_.omf) os << pop_.male_floater().size() << ' ';
    return os;
  }


  std::ostream& Simulation::stream_alleles(std::ostream& os)
  {
    collect_alleles_visitor cav;
    pop_.visit_breeder(std::ref(cav));
    os << "allele0[[length(allele0)+1]] = matrix(c("; 
    for (auto const& i : cav.v_) 
    { 
      for (auto v : i[0]) { os << v << ','; } 
    }; 
    os.seekp(-1, std::ios_base::cur) << "), nrow=" << Loci::MAX_ALLELE << ")\n";
    os << "allele1[[length(allele1)+1]] = matrix(c("; 
    for (auto const& i : cav.v_) 
    { 
      for (auto v : i[1]) { os << v << ','; } 
    }; 
    os.seekp(-1, std::ios_base::cur) << "), nrow=" << Loci::MAX_ALLELE << ")\n";
    return os;
  }


  std::ostream& Simulation::stream_mranks(std::ostream& os)
  {
    collect_mrank_visitor cbv;
    pop_.visit_breeder(std::ref(cbv));
    os << "mrank[[length(mrank)+1]] = c("; 
    for (auto const& r : cbv.v_) 
    { 
      os << r << ',';
    }; 
    os.seekp(-1, std::ios_base::cur) << ")\n";
    return os;
  }


  std::ostream& Simulation::stream_xynR(std::ostream& os)
  {
    collect_xynR_visitor cav;
    pop_.visit_patches(std::ref(cav));
    os << "xynR[[length(xynR)+1]] = matrix(c("; 
    for (auto const& x : cav.v_) 
    { 
      os << x.x << ',' << x.y << ',' << x.n << ',' << x.R << ',';
    }
    os.seekp(-1, std::ios_base::cur) << "), nrow=4)\n";
    return os;
  }


  std::ostream& Simulation::stream_takeover_stats(std::ostream& os)
  {
    os << "takeover[[length(takeover)+1]] = c(";
    auto t = (takeover_stats_ - takeover_stats_log_) / param_.log;
    os << t.attempt << ',' << t.takeover << ',' << t.walkin << ")\n";
    return os;
  }


  void run_dispatch(Parameter& param)
  {
    auto const rep = param.rep + param.repOfs;
    auto const ext = param.offile.extension();
    auto rof = param.offile;
    rof.replace_extension();
    for (size_t r = param.repOfs; r < rep; ++r)
    {
      if (rep > 1) 
      {  // adjust filename for this repetition
        auto ofn = rof;
        ofn += "_"; ofn += std::to_string(r + 1); ofn += ext;
        param.offile = ofn;
      }
      param.rep = r;
      Simulation sim(param);
      param.mode == Mating::MATING_RANDOM
        ? (param.oplacement == oPlacement::OPLACEMENT_BACK ? sim.run<Mating::MATING_RANDOM, oPlacement::OPLACEMENT_BACK>()
           : sim.run<Mating::MATING_RANDOM, oPlacement::OPLACEMENT_SORT>())
        : (param.oplacement == oPlacement::OPLACEMENT_BACK ? sim.run<Mating::MATING_RESIDENCY, oPlacement::OPLACEMENT_BACK>()
           : sim.run<Mating::MATING_RESIDENCY, oPlacement::OPLACEMENT_SORT>());
      if (param.R)
      {
        auto cmd = std::string("start ") + param.Rs + std::string(" RScript \"") + fs::absolute(param.offile).generic_string() + "\"";
        if (param.oany) std::cout << "Executing: " << cmd << '\n';
        auto err = std::system(cmd.c_str());
      }
      if (param.oany)
      {
        std::cout << "Repetition " << r + 1 << " done.\n\n"; 
      }
    }
  }

  void Run(Parameter& param)
  {
    Patch::SetupVotingSystem(param);
    run_dispatch(param);
  }

}
