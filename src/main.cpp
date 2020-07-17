/*! \file main.cpp 
 * \brief Definition of the entry point of the program
*/

#include <iostream>
#include "cmd_line.h"
#include "npm.h"      // include our model stuff


const char* NPMHelp = R"(Usage: npm [OPTION]... [OPTIONAL PARAMETER]... PARAMETER...
Options:
  --help           prints this text and exits
  -v, --verbose    prints all of the following
  -ot              prints time
  -og              prints average group size (females)
  -om              prints average group size (males)
  -off             prints number of female floater
  -omf             prints number of male floater
  -oa              prints average A0,A1,B0,B1
  -oxy             prints average x y first offspring
  -oto             prints number of takeover attempts, takeovers and walk-ins
  -aloglast        log alleles only for the last time-step

Optional parameter as name=value pairs (in brackets the default values):
  m           number of patches (1000)
  m0          percentage of initially occupied patches (90)
  nmf         initial number of male floaters (0)
  F0          baseline fecundity (1.0)
  phi         scramble competition parameter in F(n,R) (0.1)
  delta       contest competition parameter in F(n,R) (0.0)
  k           helping parameter in F(n,R) (10)
  Alleles     initial alleles A0,A1,A2,B0,B1,B2 as string ('5 0 0 5 0 0')
  Mask        alleles mask as string ('1 1 1 1 1 1')
  Sb          survival probability breeder (0.8)
  Sm          survival probability male (0.8)
  Sff         survival probability female floater (0.6)
  Smf         survival probability male floater (0.8)
  Smax        maximum survival (longevity) (0.95)
  sigma       shape parameter in Sx() (1.0)
  eps         patch search efficiency (0.005)
  t0          baseline takeover probability (0.05)
  tau         benefit for communal territory defense (1.0)
  mu          mutation probability (0.01)
  gamma       scaling parameter mutation distribution (0.1)
  ovote       'ignore' or 'account' ('account')
  bvote       one of 'ignore', 'kin', 'despotic', 'egalitarian' 'hierarchical' ('despotic')
  oplacement  offspring hierarchies placement 'back' or 'sort' ('sort')
  rep         repetitions (1)
  repOfs      start of repetition counter (0)
  R           invoke R-server with result file (false)
  Rs          start command options ('/B')
  ticks       time ticks to run (1000)
  clog        console log interval (1000)
  precision   precision of allele output (3)

Required parameter as name=value pairs:
  mode      mating mode, 'random' or 'residency'
  log       log interval, if 0 only the last state is logged
  file      file name of the result file

Examples:
  npm --verbose mode=random nmf=900 log=100 file=res1.R
  npm -v mode=residency nmf=0 eps=0.0001 mu=0.001 gamma=0.01 ticks=1e6 log=10000 file=res2.R
)";


//! The (mandatory) function that is called from the OS
//!
//! \param argc number of command line arguments
//! \param argv array of command line arguments (C-strings)
int main(int argc, const char* argv[])
{
  // Get command line arguments
  cmd::cmd_line_parser clp(argc, argv);

  // just --help ? 
  if (clp.flag("--help"))
  {
    std::cout << NPMHelp;
    return 0;
  }

  //
  // Fine, someone want to run a simulation
  //

  try
  {
    npm::Parameter param;                         // default parameter set

    // Console output flags
    param.verbose = clp.flag("-v") || clp.flag("--verbose");
    param.ot = clp.flag("-ot") || param.verbose;
    param.og = clp.flag("-og") || param.verbose;
    param.om = clp.flag("-om") || param.verbose;
    param.off = clp.flag("-off") || param.verbose;
    param.omf = clp.flag("-omf") || param.verbose;
    param.oa = clp.flag("-oa") || param.verbose;
    param.oxy = clp.flag("-oxy") || param.verbose;
    param.oto = clp.flag("-oto") || param.verbose;
    param.aloglast = clp.flag("-aloglast") || param.aloglast;
    param.oprof = clp.flag("-prof");
    param.oany = param.ot || param.og || param.om || param.off || param.omf || param.oa || param.oxy || param.oto || param.oprof;

    param.log = clp.required<size_t>("log");
    param.offile = clp.required<fs::path>("file");
    if (param.offile.parent_path().empty()) {
      param.offile = fs::path(".") / param.offile;
    }
    param.nmf = clp.required<size_t>("nmf");

    std::string pstr = clp.required<std::string>("mode");
    param.mode = (npm::Mating) cmd::check_any(pstr, npm::mating_name, "invalid mode parameter");

    pstr = npm::ovote_name[(int)param.ovote];
    clp.optional("ovote", pstr);
    param.ovote = (npm::oVote)cmd::check_any(pstr, npm::ovote_name, "invalid ovote parameter");

    pstr = npm::bvote_name[(int)param.bvote];
    clp.optional("bvote", pstr);
    param.bvote = (npm::bVote)cmd::check_any(pstr, npm::bvote_name, "invalid bvote parameter");

    pstr = npm::oplacement_name[(int)param.oplacement];
    clp.optional("oplacement", pstr);
    param.oplacement = (npm::oPlacement)cmd::check_any(pstr, npm::oplacement_name, "invalid oplacement parameter");

    clp.optional("m", param.m);
    clp.optional("m0",param.m0);
    clp.optional("F0", param.F0);
    clp.optional("phi", param.phi);
    clp.optional("delta", param.delta);
    clp.optional("k", param.k);
    clp.optional("Alleles", param.alleles);
    clp.optional("Mask", param.mask);
    clp.optional("Sb", param.Sb);
    clp.optional("Sm", param.Sm);
    clp.optional("Sff", param.Sff);
    clp.optional("Smf", param.Smf);
    clp.optional("Smax", param.Smax);
    clp.optional("sigma", param.sigma);
    clp.optional("eps", param.eps);
    clp.optional("t0", param.t0);
    clp.optional("tau", param.tau);
    clp.optional("mu", param.mu);
    clp.optional("gamma", param.gamma);
    clp.optional("rep", param.rep);
    clp.optional("repOfs", param.repOfs);
    param.R = clp.flag("-R");
    clp.optional("Rs", param.Rs);
    // double? - allow 10^5 etc...
    double ticks = static_cast<double>(param.ticks);
    clp.optional("ticks", ticks);
    param.ticks = static_cast<size_t>(ticks);
    clp.optional("clog", param.clog);
    clp.optional("precision", param.precision);
    clp.optional<size_t>("nmf", param.nmf);
    // finally we can start the model...
    npm::Run(param);
    std::cout << "Regards\n";
    return 0;
  }
  catch (cmd::parse_error& e)
  {
    std::cerr << "npm: Invalid arguments: " << e.what() << '\n';
    std::cerr << "use 'npm --help'\nfor instructions.\n";
  }
  catch (std::exception& e)
  {
    std::cerr << "npm: Fatal error: " << e.what() << '\n';
  }
  catch (...)
  {
    std::cerr << "npm: Fatal error: unknown exception\n";
  }
  return -1;     // return error code
}
