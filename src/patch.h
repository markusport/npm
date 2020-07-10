/*! \file patch.h
* \brief Declaration of a patch (home-range) in the natal philopatry model
*
*/

#ifndef NPM_PATCH_H_INCLUDED
#define NPM_PATCH_H_INCLUDED

#include <vector>
#include <cmath>
#include <memory>
#include "individual.h"


namespace npm {


  //! \brief behavior record of an offspring
  struct xynR_type
  {
    double x;   //!< voted x
    double y;   //!< voted y
    size_t n;   //!< group size
    size_t R;   //!< mothers rank
  };


  // container type for individuals
  typedef std::vector<Individual> container_t;


  //! \brief A patch
  //!
  //! A patch is a collection of breeders. A.k.a home-range
  //! The first breeder in the collection is the dominant one.
  class Patch
  {
  public:
    static void SetupVotingSystem(Parameter const& param);

    //! \brief creates an empty patch
    Patch() = default;

    //! \brief creates an occupied patch
    //! \param dominant The individual that becomes dominant
    //! \param male The individual that becomes male breeder (could be 0)
    explicit Patch(Individual const& dominant, const Individual* male);
    
    //! \brief Returns true if the patch is empty
    bool empty() const { return breeder_.empty(); }

    //! Returns number of breeders
    size_t size() const { return breeder_.size(); }

    //! \brief Returns an pointer to the breeding male
    //! If the returned pointer is nullptr, no breeding male
    //! exist on the patch.
    Individual const* male() const { return male_.empty() ? nullptr : &*male_.cbegin(); }

    //! \brief Returns an pointer to the breeding male
    //! If the returned pointer is nullptr, no breeding male
    //! exist on the patch.
    Individual* male() { return male_.empty() ? nullptr : &*male_.begin(); }

    //! \brief Sets new male
    //! \param newMale The new male
    void set_male(Individual const& newMale) { male_.assign(1, newMale); }

    //! \brief Returns the breeder collection. 
    container_t const& breeder() const { return breeder_; }

    //! \brief Returns an iterator to the first breeder_
    //! If the returned iterator is equal to cend(), the patch
    //! is empty.
    container_t& breeder() { return breeder_; }

    //! \brief Returns the rank vector of the female offspring
    std::vector<xynR_type> const& verdict() const { return verdict_; }

    //! brief Handles survival on the patch
    //! \param param parameter set
    void do_survival(Parameter const& param);

    //! \brief Handles reproduction on patch
    //! \tparam MODE Mode
    //! \param param parameter set
    //! \param male_floater male floater pool
    template <Mating MODE>
    void do_reproduction(Parameter const& param, container_t const& male_floater);

    //! \brief Handles dispersal on patch and to the floater pool
    //! \tparam PLACEMENT oPlacemanet
    //! \param param parameter set
    //! \param female_floater the female floater pool
    //! \param male_floater the male floater pool
    template <oPlacement PLACEMENT>
    void do_dispersal(Parameter const& param, container_t& female_floater, container_t& male_floater);

    //! \brief Handles colonization of the patch by female floater
    //! \param param parameter set
    //! \param floater The female floater
    void do_colonization(Parameter const& param, Individual const& floater);

  private:
    void prepare_reproduction();
    void create_offsprings(Parameter const& param, Individual const& male);
    void disperse_males_and_poll(container_t& male_floater);

    // voting system
    template <oVote V> double offspring_vote(size_t ioffs) const;
    template <bVote V> double breeder_vote(size_t ioffs) const;

    container_t breeder_;
    container_t female_offspring_;
    container_t male_offspring_;
    container_t male_;
    std::vector<double> x_;             // offspring stay prob
    std::vector<double> y_;             // mothers accept prob
    std::vector<unsigned> R_;           // rank of mother == position of mother + 1
    std::vector<xynR_type> verdict_;    // the outcome of the poll x,y,n,R
  };

}

#endif
