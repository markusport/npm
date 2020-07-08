/*! \file population.h
* \brief Declaration of the population in the natal philopatry model
*
*/

#ifndef NPM_POPULATION_H_INCLUDED
#define NPM_POPULATION_H_INCLUDED

#include "patch.h"


namespace npm{


  //! \brief The population
  //!
  //! A Population consists of its patches and the floater pool.
  class Population
  {
  public:
    //! creates a empty population
    Population() = default;

    //! \brief creates the initial population
    //! \param param parameter set
    //!
    //! Creates param.m patches, param.m0% of them occupied
    //! with a breeder pair.
    explicit Population(Parameter const& param);

    //! Returns the patches
    std::vector<Patch> const& patches() const { return patches_; }

    //! Returns the patches, non const
    std::vector<Patch>& patches() { return patches_; }

    //! Returns the female floaters
    container_t const& female_floater() const { return female_floater_; }

    //! Returns the female floaters, non const
    container_t& female_floater() { return female_floater_; }

    //! Returns the male floaters
    container_t const& male_floater() const { return male_floater_; }

    //! Returns the male floaters, non const
    container_t& male_floater() { return male_floater_; }

    //! \brief Random shuffle of floaters
    void shuffle_floater(Parameter const& param);

    //! \brief Handles survival of the floater
    //! \param param parameter set
    void do_floater_survival(Parameter const& param);

    //! \brief Handles colonization and takeover
    //! \tparam MODE Mode::RANDOM_MATING or Mode::MALE_RESIDENCY
    //! \param param parameter set
    //! \returns { number of takeover attempts, number of takeovers }
    template <Mating MODE>
    TakeoverStats do_colonization(Parameter const& param);

    //! \brief
    //! \param fun a function object with the signature void fun(Individual const&);
    //!
    //! Applies the given function object to every individual in
    //! this population.
    template <typename UnaryFunction>
    void visit_all(UnaryFunction fun);


    //! \brief
    //! \param fun a function object with the signature void fun(Individual const&);
    //!
    //! Applies the given function object to every breeder in
    //! this population.
    template <typename UnaryFunction>
    void visit_breeder(UnaryFunction fun);


    //! \brief
    //! \param fun a function object with the signature void fun(Individual const&);
    //!
    //! Applies the given function object to every patch in
    //! this population.
    template <typename UnaryFunction>
    void visit_patches(UnaryFunction fun);


  private:
    std::vector<Patch> patches_;
    container_t female_floater_;
    container_t male_floater_;
  };

  

  //
  // declaration of template specializations
  //

  template <>
  TakeoverStats Population::do_colonization<Mating::RANDOM>(Parameter const& param);


  template <>
  TakeoverStats Population::do_colonization<Mating::RESIDENCY>(Parameter const& param);


  //
  // implementation of template member functions
  //

  template <typename UnaryFunction>
  inline void Population::visit_all(UnaryFunction fun)
  {
    for (auto& patch : patches_)
    {
      for (auto& ind : patch.breeder()) fun(ind);
      if (patch.male()) fun(*patch.male());
    }
    for(auto& ind : female_floater_) fun(ind);
    for(auto& ind : male_floater_) fun(ind);
  }


  template <typename UnaryFunction>
  inline void Population::visit_breeder(UnaryFunction fun)
  {
    for (auto& patch : patches_)
    {
      for (auto& ind : patch.breeder()) fun(ind);
    }
  }


  template <typename UnaryFunction>
  inline void Population::visit_patches(UnaryFunction fun)
  {
    for (auto const& patch : patches_)
    {
      if (!patch.empty())
      {
        fun(patch);
      }
    }
  }


}

#endif
