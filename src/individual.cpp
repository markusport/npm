/*! \file individual.cpp
* \brief Definition of an individual in the natal philopatry model
*/

#include "rndutils.hpp"
#include "individual.h"


namespace npm {


  Individual::Individual(Parameter const& param)
  {
    phen = inherited[0] = inherited[1] = param.alleles;
    age = 0;
    mRank = 0;
  }


  Individual::Individual(Parameter const& param, Individual const& female, Individual const& male, unsigned mRank)
  {
    age = 0;
    this->mRank = mRank;

    std::bernoulli_distribution bernoulli_mu(param.mu);
    mutation_dist rndMut(0.0, param.sigma);

    // Recombination - 8 bit of randomness
    rndutils::binary_distribution binary_dist;
    for (size_t i = 0; i < Loci::MAX_ALLELE; ++i)
    {
      // Recombination
      auto recomb = binary_dist(RndEng);
      auto x = female.inherited[recomb][i];
      auto y = male.inherited[recomb][i];
      // Mutation
      if (bernoulli_mu(RndEng)) x += rndMut(RndEng);
      if (bernoulli_mu(RndEng)) y += rndMut(RndEng);
      // optional: mask out unused alleles
      inherited[0][i] = param.mask[i] * x;
      inherited[1][i] = param.mask[i] * y;
      phen[i] = 0.5 * (x + y);
    }
  }


}
