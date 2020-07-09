/*! \file visitors.h
* \brief A growing collection of population visitors
*
*/

#ifndef NPM_VISITORS_H_INCLUDED
#define NPM_VISITORS_H_INCLUDED

#include <ostream>
#include <utility>
#include <functional>
#include "individual.h"
#include "patch.h"


namespace npm{

  struct xynR_collection
  {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<size_t> n;
    std::vector<size_t> R;
  };


  //! \brief collect x y n R visitor
  class collect_xynR_visitor
  {
  public:
    collect_xynR_visitor() {}

    void operator()(Patch const& patch)
    {
      v_.insert(v_.end(), patch.verdict().begin(), patch.verdict().end());
    }
    std::vector<xynR_type> v_;
  };


  //! \brief collect mothers rank
  class collect_mrank_visitor
  {
  public:
    collect_mrank_visitor() {}

    void operator()(Individual const& x)
    {
      v_.push_back(x.mRank);
    }
    std::vector<unsigned> v_;
  };


  //! \brief collect alleles visitor
  class collect_alleles_visitor
  {
  public:
    collect_alleles_visitor() {}

    void operator()(Individual const& x)
    {
      v_.push_back(x.inherited);
    }

    std::vector<std::array<Alleles, 2>> v_;
  };


  //! \brief mean allele visitor
  class mean_allele_visitor
  {
  public:
    mean_allele_visitor() : counts_(0) 
    {
      alleles_ = { 0.0 };
    }
    
    void operator()(Individual const& x)
    {
      if (counts_)
      {
        for (size_t i=0; i < Loci::MAX_ALLELE; ++i) alleles_[i] += x.phen[i];
      }
      else
      {
        for (size_t i = 0; i < Loci::MAX_ALLELE; ++i) alleles_[i] = x.phen[i];
      }
      ++counts_;
    }
    
    Alleles mean() const
    {
      if (counts_ == 0) return alleles_;
      Alleles tmp(alleles_);
      for (size_t i = 0; i < Loci::MAX_ALLELE; ++i) tmp[i] /= counts_;
      return tmp;
    }

  private:
    size_t counts_;
    Alleles alleles_;
  };


  inline xynR_type& operator+=(xynR_type& lhs, xynR_type const& rhs)
  {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.n += rhs.n;
    lhs.R += rhs.R;
    return lhs;
  }


  //! \brief mean behavior first offspring visitor
  class mean_behavior_visitor
  {
  public:
    mean_behavior_visitor() : c_(0), sum_{0,0,0,0}
    {
    }

    void operator()(Patch const& patch)
    {
      if (!patch.verdict().empty()) 
      { 
        auto xynR = patch.verdict()[0]; 
        if (c_ == 0) 
        {
          sum_ = xynR; 
        }
        else 
        {
          sum_ += xynR;
        }
        ++c_;
      }
    }

    xynR_type mean() const
    {
      if (c_ == 0) return {0,0,0,0};
      return { sum_.x / c_, sum_.y / c_, sum_.n / c_, sum_.R / c_ };
    }

  private:
    xynR_type sum_;
    size_t c_;
  };

}

#endif
