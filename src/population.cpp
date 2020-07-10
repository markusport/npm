/*! \file population.cpp
* \brief Definition the population in the natal philopatry model
*/

#include <algorithm>
#include "population.h"


namespace npm {


  Population::Population(Parameter const& param)
  {
    auto M = static_cast<size_t>(std::ceil((param.m0 * param.m / 100.0)));
    auto Default = Individual(param);
    for (size_t i = 0; i < M; ++i) 
    {
      patches_.emplace_back(Default, (param.mode==Mating::MATING_RESIDENCY ? &Default : nullptr));
    }
    for (size_t i = M; i < param.m; ++i) patches_.emplace_back();
    for (size_t i=0; i < param.nmf; ++i) male_floater_.emplace_back(Default);
  }


  void Population::shuffle_floater(Parameter const& param)
  {
    std::shuffle(female_floater_.begin(), female_floater_.end(), RndEng);
    std::shuffle(male_floater_.begin(), male_floater_.end(), RndEng);
  }


  void Population::do_floater_survival(Parameter const& param)
  {
    std::bernoulli_distribution bernoulli_Sff(1.0 - param.Sff);
    female_floater_.erase(std::remove_if(female_floater_.begin(), female_floater_.end(), [&bernoulli_Sff](Individual&)
    {
      return bernoulli_Sff(RndEng);
    }), female_floater_.end());
    std::bernoulli_distribution bernoulli_Smf(1.0 - param.Smf);
    male_floater_.erase(std::remove_if(male_floater_.begin(),male_floater_.end(),[&bernoulli_Smf](Individual&)
    {
      return bernoulli_Smf(RndEng);
    }),male_floater_.end());
  }


  template <>
  TakeoverStats Population::do_colonization<Mating::MATING_RANDOM>(Parameter const& param)
  {
    TakeoverStats tc{0, 0, 0};
    if (female_floater_.empty()) return tc;
    std::poisson_distribution<> rndPois(param.eps * female_floater_.size());
    for (auto& patch : patches_)
    {
      if (female_floater_.empty()) return tc;
      int k = rndPois(RndEng);
      tc.attempt += k;
      if (k)
      {
        bool takeover = false;
        if (patch.empty())
        {
          ++tc.walkin;
          takeover = true;
        }
        else
        {
          double tprob = static_cast<double>(k) * param.t0 * std::exp(-param.tau * (patch.size() - 1));
          takeover = std::bernoulli_distribution(tprob)(RndEng);
        }
        if (takeover)
        {
          patch.do_colonization(param, female_floater_.back());
          female_floater_.pop_back();
          ++tc.takeover;
        }
      }
    }
    return tc;
  }


  template <>
  TakeoverStats Population::do_colonization<Mating::MATING_RESIDENCY>(Parameter const& param)
  {
    auto tc = do_colonization<Mating::MATING_RANDOM>(param);  // same for females
    for (auto& patch : patches_)
    {
      if (male_floater_.empty()) return tc;
      if (nullptr == patch.male())
      {
        patch.set_male(male_floater_.back());
        male_floater_.pop_back();
      }
    }
    return tc;
  }


}