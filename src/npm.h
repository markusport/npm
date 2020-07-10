/*! \file npm.h
 * \brief Public declarations of the natal philopatry model
 *
*/

#ifndef NPM_NPM_H_INCLUDED
#define NPM_NPM_H_INCLUDED

#include <array>
#include <string>
#include <filesystem>
#include "rndutils.hpp"


namespace fs = std::filesystem;


//! \brief Everything of the Natal philopatry model
namespace npm {


  extern const char* Version;


  //! \brief random number engine
  extern rndutils::xorshift128 thread_local RndEng;


  //! \brief Mutation distribution
  using mutation_dist = std::cauchy_distribution<>;


  //! \brief mating modes
  enum Mating
  {
    MATING_RANDOM,
    MATING_RESIDENCY,
    MATING_MAX
  };


  //! \brief offspring vote
  enum oVote
  {
    OVOTE_IGNORE,         //!< always 1.0
    OVOTE_ACCOUNT,        //!< x(n,R)
    OVOTE_MAX
  };


  //! \brief breeder vote
  enum bVote
  {
    BVOTE_IGNORE,         //!< always 1.0
    BVOTE_KIN,            //!< y(n,R)
    BVOTE_DESPOTIC,       //!< y(n,1)
    BVOTE_EGALITARIAN,    //!< average y(n,R)
    BVOTE_HIERARCHICAL,   //!< average y(n,R), R <= R mother
    BVOTE_MAX
  };


  //! \brief offspring placement
  enum oPlacement
  {
    OPLACEMENT_BACK,
    OPLACEMENT_SORT,
    OPLACEMENT_MAX
  };


  extern const char* mating_name[Mating::MATING_MAX];
  extern const char* ovote_name[oVote::OVOTE_MAX];
  extern const char* bvote_name[bVote::BVOTE_MAX];
  extern const char* oplacement_name[oPlacement::OPLACEMENT_MAX];
  

  //! \brief allele gene loci
  enum Loci
  {
    A0, A1, A2, B0, B1, B2, MAX_ALLELE
  };


  //! \brief A set of alleles, e.g. {A0,A1,A2,...,B0,B1,B2,...}
  using Alleles = std::array<double, MAX_ALLELE>;


  //! \brief Helper class to track takeovers
  struct TakeoverStats
  {
    size_t attempt;   //! < counter of attempts
    size_t takeover;  //! < counter of successful takeovers
    size_t walkin;    //! < counter of colonization of empty patches
  
    TakeoverStats& operator+=(TakeoverStats const& rhs)
    {
      attempt += rhs.attempt;
      takeover += rhs.takeover;
      walkin += rhs.walkin;
      return *this;
    }

    TakeoverStats operator-(TakeoverStats const& rhs) const
    {
      return { attempt - rhs.attempt, takeover - rhs.takeover, walkin - rhs.walkin };
    }

    TakeoverStats operator/(size_t rhs) const
    {
      return { attempt / rhs, takeover / rhs, walkin / rhs };
    }
  };


  //! \brief Parameter set of the model
  struct Parameter
  {
    size_t m = 1000;                          //!< number of patches
    double m0 = 90;                           //!< initially occupied patches [%]
    size_t nmf = 0;                           //!< initial number of male floaters
    size_t F0 = 1;                            //!< baseline fecundity
    double phi = 0.1;                       //!< parameter in F(n,R)
    double delta = 0.0;                       //!< parameter in F(n,R)
    double k = 10.0;                          //!< parameter in F(n,R)
    Alleles alleles = { 
      5.0, 0.0, 0.0,  5.0, 0.0, 0.0 
    };                                        //!< initial alleles
    Alleles mask = {
      1.0, 1.0, 1.0,  1.0, 1.0, 1.0
    };                                        //!< masking factor
    double Sb = 0.8;                          //!< baseline survival probability breeder
    double Sm = 0.8;                          //!< baseline survival probability male
    double Sff = 0.6;                         //!< survival probability female floater
    double Smf = 0.8;                         //!< survival probability male floater
    double Smax = 0.95;                       //!< maximum survival (longevity)
    double sigma = 1.0;                       //!< parameter in Sx(n)
    double eps = 0.005;                       //!< patch search efficiency
    double t0 = 0.05;                         //!< baseline takeover probability
    double tau = 1.0;                           //!< benefit for communal territory defense
    double mu = 0.1;                          //!< mutation probability
    double gamma = 0.01;                      //!< standard deviation mutation distribution
    size_t ticks = 1000;                      //!< Time ticks to run
    size_t rep = 1;                           //!< Repetitions
    size_t repOfs = 0;                        //!< Start of repetition counter
    bool R = false;                           //!< invoke R server with result file
    std::string Rs = "/B";                    //!< R start command
    size_t log = 0;                           //!< log interval
    size_t clog = 1000;                       //!< console log interval
    bool aloglast = false;                    //!< if true, log alleles for last timestep only
    unsigned precision = 3;                   //!< precision of allele output
    fs::path offile = "";                     //!< output data file
    bool verbose = false;                     //!< verbose output
    bool ot = false;                          //!< print time 
    bool og = false;                          //!< print average group size
    bool om = false;                          //!< print average number of males
    bool off = false;                         //!< print number of female floaters
    bool omf = false;                         //!< print number of male floaters
    bool oa = false;                          //!< print average {A0,A1,B0,B1}
    bool oxy = false;                         //!< print average x y first offspring
    bool oto = false;                         //!< print takeover stats
    bool oprof = false;                       //!< profiling
    bool oany = false;                        //!< any of the above
    Mating mode;                              //!< mating mode
    oVote ovote = oVote::OVOTE_ACCOUNT;             //!< offspring vote
    bVote bvote = bVote::BVOTE_DESPOTIC;            //!< breeder vote
    oPlacement oplacement = oPlacement::OPLACEMENT_SORT; //!< offspring placement mode


    double thetaB() const { return (Sb - Smax * (1.0 - std::exp(-sigma))) / std::exp(-sigma); }
    double thetaM() const { return (Sm - Smax * (1.0 - std::exp(-sigma))) / std::exp(-sigma); }
  };


  //! \brief Run the model
  //! \param param the parameter set
  void Run(Parameter& param);


}


#endif
