# npm
Natal pylopatry model

```
Usage: npm [OPTION]... [OPTIONAL PARAMETER]... PARAMETER...
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
  alpha       parameter in F(n,R) (0.1)
  delta       parameter in F(n,R) (0.0)
  k           parameter in F(n,R) (10)
  Alleles     initial alleles as string ('5 0 0 5 0 0')
  Mask        alleles mask as string ('1 1 1 1 1 1')
  Sb          survival probability breeder (0.8)
  Sm          survival probability male (0.8)
  Sff         survival probability female floater (0.6)
  Smf         survival probability male floater (0.8)
  Smax        maximum survival (longevity) (0.95)
  gamma       shape parameter in Sx() (1.0)
  eps         patch search efficiency (0.005)
  t0          baseline takeover probability (0.05)
  d           benefit for communal territory defense (1.0)
  mu          mutation probability (0.01)
  sigma       scaling parameter mutation distribution (0.1)
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
  npm -v mode=residency nmf=0 eps=0.0001 mu=0.001 sigma=0.01 ticks=1e6 log=10000 file=res2.R

```
