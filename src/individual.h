/*! \file individual.h
* \brief Declaration of an individual in the natal philopatry model.
*
*/

#ifndef NPM_INDIVIDUAL_H_INCLUDED
#define NPM_INDIVIDUAL_H_INCLUDED

#include "npm.h"


namespace npm {


  //! \brief An individual
  //!
  //! An Individual is not much more than a bag of its alleles
  struct Individual
  {
    // \brief default copy constructor and assignment operator are fine
    Individual(Individual const&) = default;
    Individual& operator=(Individual const&) = default;
    
    //! \brief creates individual with undefined state
    Individual() {};
    
    //! \brief creates individual with initial alleles
    explicit Individual(Parameter const& param);

    //! \brief Creates an offspring
    //! \param param the parameter set
    //! \param female Female ancestor
    //! \param male male ancestor
    //! \param mRank mothers rank
    //!
    //! Mutation and recombination happens here
    Individual(Parameter const& param, Individual const& female, Individual const& male, unsigned mRank);

    Alleles phen;						            //!< active 'phenotype'
    std::array<Alleles, 2> inherited;   //!< inherited alleles [mother, father]
    unsigned age;                       //!< age of this individual
    unsigned mRank;                     //!< mothers rank
  };


}

#endif
