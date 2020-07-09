/*! \file patch.cpp
* \brief Definition of a patch in the natal philopatry model
*/

#include <random>
#include <algorithm>
#include "patch.h"


namespace npm {

  namespace {

    // member function pointers of the voting system
    double (Patch::*offspring_vote_pmf)(size_t) const;
    double (Patch::*breeder_vote_pmf)(size_t) const;

  }


  //! \brief patch fecundity probability
  struct fecundity
  {
    //! \param param the parameter set
    //! \param n number of occupied sites
    fecundity(Parameter const& param, double n)
    : fact( param.F0 * (1.0 - param.alpha * n) * (1.0 - std::exp(-param.k * n)) ),
      delta(param.delta)
    {}

    //! \param Rank rank of the mother in spe
    //! \returns F(n,R)
    double operator()(double Rank) const 
    {
      return fact * std::pow(Rank, -delta);
    }

    const double fact;
    const double delta;
  };


  //! \brief Probability for offspring to stay on a patch
  //! \param x an individual
  //! \param n number of occupied sites
  //! \param R rank of the mother
  //!
  //! \returns x(n, R)
  inline double stay_probability(Individual const& x, double n, double R)
  {
    return 1.0 / (1.0 + std::exp(x.phen[B0] + n * x.phen[B1] + R * x.phen[B2]));
  }


  //! \brief Probability to accept an offspring
  //! \param x an individual
  //! \param n number of occupied sites
  //! \param R rank of the mother
  //!
  //! \returns y(n, R)
  inline double accept_probability(Individual const& x, double n, double R)
  {
    return 1.0 / (1.0 + std::exp(x.phen[A0] + n * x.phen[A1] + R * x.phen[A2]));
  }


  Patch::Patch(Individual const& dominant, Individual const* male)
  {
    breeder_.assign(1, dominant);
    if (male) set_male(*male);
  }


  void do_mortality(container_t& c, double SurvivalProp)
  {
    std::bernoulli_distribution bernoulli_P(1.0 - SurvivalProp);
    c.erase(std::remove_if(c.begin(), c.end(), [&bernoulli_P](Individual& ind)
    {
      return (ind.age > 0) && (bernoulli_P(RndEng));
    }), c.end());
  }


  void Patch::do_survival(Parameter const& param)
  {
    auto n = static_cast<double>(breeder_.size());
    auto thetaB = param.thetaB(); 
    auto thetaM = param.thetaM(); 
    do_mortality(breeder_, thetaB + (param.Smax - thetaB) * (1.0 - std::exp(-param.gamma * n)));
    do_mortality(male_, thetaM + (param.Smax - thetaM) * (1.0 - std::exp(-param.gamma * n)));
  }


  template <>
  void Patch::do_reproduction<Mating::MATING_RANDOM>(Parameter const& param, container_t const& male_floater)
  {
    prepare_reproduction();
    if (!(male_floater.empty() || empty()))
    {
      // select male at random
      rndutils::uniform_signed_distribution<> rndMale(0, (int)male_floater.size() - 1);
      Individual const* male = &male_floater[rndMale(RndEng)];
      create_offsprings(param, *male);
    }
  }
  

  template <>
  void Patch::do_reproduction<Mating::MATING_RESIDENCY>(Parameter const& param,container_t const&)
  {
    prepare_reproduction();
    if (!(male_.empty() || empty()))
    {
      create_offsprings(param, male_[0]);
    }
  }


  void Patch::prepare_reproduction()
  {
    female_offspring_.clear();
    male_offspring_.clear();
    x_.clear();
    y_.clear();
    R_.clear();
    verdict_.clear();
  }


  void Patch::create_offsprings(Parameter const& param, Individual const& male)
  {
    rndutils::binary_distribution binary_dist;
    double n = static_cast<double>(breeder_.size());
    auto fecundityP = fecundity(param, n);
    for (size_t i = 0; i < breeder_.size(); ++i)
    {
      auto R = static_cast<double>(i + 1);
      y_.push_back( accept_probability(breeder_[i], n, R) );
      std::bernoulli_distribution bernoulli_fecundity( fecundityP(R) );
      for (int k=0; k<param.F0; ++k) 
      {
        if (bernoulli_fecundity(RndEng))
        {
          if (binary_dist(RndEng))
          { // female offspring
            female_offspring_.emplace_back(param, breeder_[i], male, static_cast<unsigned>(i + 1));
            x_.push_back(stay_probability(female_offspring_.back(), n, R));
            R_.push_back(static_cast<unsigned>(i + 1));
          }
          else
          {  // male offspring
            male_offspring_.emplace_back(param, breeder_[i], male, static_cast<unsigned>(i + 1));
          }
        }
      }
    }
  }


  void Patch::disperse_males_and_poll(container_t& male_floater)
  {
    // males goes en block to the male floater pool.
    male_floater.insert(male_floater.end(), male_offspring_.begin(), male_offspring_.end());
    if (empty()) 
    { // lonely males become floaters again
      male_floater.insert(male_floater.end(), male_.begin(), male_.end());
      male_.clear();
      return;
    }

    // perform the poll
    auto const oldN = breeder_.size();
    for (size_t i = 0; i < female_offspring_.size(); ++i)
    {
      verdict_.push_back({ (this->*offspring_vote_pmf)(i), (this->*breeder_vote_pmf)(i), oldN, R_[i] });
    }
  }


  template <>
  void Patch::do_dispersal<oPlacement::OPLACEMENT_BACK>(Parameter const& param, container_t& female_floater, container_t& male_floater)
  {
    disperse_males_and_poll(male_floater);
    auto const oldN = breeder_.size();
    for (size_t i = 0; i < female_offspring_.size(); ++i)
    {
      if (std::bernoulli_distribution(verdict_[i].x * verdict_[i].y)(RndEng))
      { // stay on patch
        breeder_.emplace_back(female_offspring_[i]);
      }
      else
      { // emigrate to floater pool
        female_floater.push_back(female_offspring_[i]);
      }
    }
    // Shuffle new breeders
    std::shuffle(breeder_.begin() + oldN, breeder_.end(), RndEng);
  }


  template <>
  void Patch::do_dispersal<oPlacement::OPLACEMENT_SORT>(Parameter const& param, container_t& female_floater, container_t& male_floater)
  {
    disperse_males_and_poll(male_floater);
    size_t rank_shift = 0;
    for (size_t i = 0; i < female_offspring_.size(); ++i)
    {
      if (std::bernoulli_distribution(verdict_[i].x * verdict_[i].y)(RndEng))
      { // stay on patch
        auto pos = breeder_.begin() + (R_[i] + rank_shift++);
        breeder_.insert(pos, female_offspring_[i]);
      }
      else
      { // emigrate to floater pool
        female_floater.push_back(female_offspring_[i]);
      }
    }
  }


  void Patch::do_colonization(Parameter const& param, Individual const& floater)
  {
    breeder_.assign(1, floater);
  }


  // Voting system
  
  template <> double Patch::offspring_vote<oVote::OVOTE_IGNORE>(size_t) const
  {
    return 1.0;
  }


  template <> double Patch::offspring_vote<oVote::OVOTE_ACCOUNT>(size_t ioffs) const
  {
    return x_[ioffs];
  }


  template <> double Patch::breeder_vote<bVote::BVOTE_IGNORE>(size_t) const
  {
    return 1.0;
  }


  template <> double Patch::breeder_vote<bVote::BVOTE_KIN>(size_t ioffs) const
  {
    return y_[R_[ioffs] - 1];
  }


  template <> double Patch::breeder_vote<bVote::BVOTE_DESPOTIC>(size_t) const
  {
    return y_[0];
  }


  template <> double Patch::breeder_vote<bVote::BVOTE_EGALITARIAN>(size_t) const
  {
    double res = 0.0;
    auto N = breeder_.size();
    for (size_t j = 0; j < N; ++j)
    {
      res += y_[j];
    }
    return res / N;
  }


  template <> double Patch::breeder_vote<bVote::BVOTE_HIERARCHICAL>(size_t ioffs) const
  {
    double res = 0.0;
    auto N = static_cast<size_t>(R_[ioffs]);
    for (size_t j = 0; j < N; ++j)
    {
      res += y_[j];
    }
    return res / N;
 }


  void Patch::SetupVotingSystem(Parameter const& param)
  {
    switch (param.ovote)
    {
      case oVote::OVOTE_IGNORE: offspring_vote_pmf = &Patch::offspring_vote<oVote::OVOTE_IGNORE>; break;
      case oVote::OVOTE_ACCOUNT: offspring_vote_pmf = &Patch::offspring_vote<oVote::OVOTE_ACCOUNT>; break;
    }
    switch (param.bvote)
    {
      case bVote::BVOTE_IGNORE: breeder_vote_pmf = &Patch::breeder_vote<bVote::BVOTE_IGNORE>; break;
      case bVote::BVOTE_KIN: breeder_vote_pmf = &Patch::breeder_vote<bVote::BVOTE_KIN>; break;
      case bVote::BVOTE_DESPOTIC: breeder_vote_pmf = &Patch::breeder_vote<bVote::BVOTE_DESPOTIC>; break;
      case bVote::BVOTE_EGALITARIAN: breeder_vote_pmf = &Patch::breeder_vote<bVote::BVOTE_EGALITARIAN>; break;
      case bVote::BVOTE_HIERARCHICAL: breeder_vote_pmf = &Patch::breeder_vote<bVote::BVOTE_HIERARCHICAL>; break;
    }
  }


}

