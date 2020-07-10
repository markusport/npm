# npm
Natal pylopatry model


## Build the simulation

If you are on Windows, you *can* skip the build-process and use the precompiled binary `npm.exe` in the `bin` folder. Otherwise, go ahead and build it with CMake:

```
:~npm$ mkdir build && cd build
:~npm/build$ cmake ..
:~npm/build$ cmake --build . --config Release --target install
```

Tested on Linux (g++ > 8.0), MacOS (Xcode > 10) and Windows (Visual Studio 2019), this should have created the binary `:~npm/bin/npm`:

```
:~npm/bin ./npm --help
Usage: npm [OPTION]... [OPTIONAL PARAMETER]... PARAMETER...
Options:
  --help           prints this text and exits
  -v, --verbose    prints all of the following
  -ot              prints time
  -og              prints average group size (females)
  -om              prints average group size (males)
  -off             prints number of female floater
  -omf             prints number of male floater
  -oa              prints average A0,A1,A2,B0,B1,B2; alpha0, alpha1, beta0, beta1 in eq. (2a,b) in Port et al.; A2 and B2 model dependence of x and y, respectively,      on mother's dominance rank. Not addressed in Port et al.
  -oxy             prints average x y first offspring
  -oto             prints number of takeover attempts, takeovers and walk-ins
  -aloglast        log alleles only for the last time-step

Optional parameter as name=value pairs (in brackets the default values):
  m           number of patches (1000); patches referred to "home-ranges" in Port et al. 
  m0          percentage of initially occupied patches (90); patches referred to "home-ranges" in Port et al. 
  nmf         initial number of male floaters (0)
  F0          baseline fecundity (1.0)
  alpha       scramble competition parameter in F(n,R) (0.1); phi in Port et al. 
  delta       contest competition parameter in F(n,R) (0.0)
  k           helping parameter in F(n,R) (10); not adressed in Port et al. 
  Alleles     initial alleles as string ('5 0 0 5 0 0')
  Mask        alleles mask as string ('1 1 1 1 1 1'); incorporates option to mask effect of A0,A1,A2,B0,B1,B2 on x and y. Mask is 1 1 0 1 1 0 in Port et al. since effect of mother's rank (A2, B2) on x and y is not modelled
  Sb          survival probability breeder (0.8)
  Sm          survival probability male (0.8)
  Sff         survival probability female floater (0.6)
  Smf         survival probability male floater (0.8)
  Smax        maximum survival (longevity) (0.95)
  gamma       shape parameter in Sx() (1.0); sigma in Port et al.
  eps         patch search efficiency (0.005)
  t0          baseline takeover probability (0.05)
  d           benefit for communal territory defense (1.0); tau in Port et al.
  mu          mutation probability (0.01)
  sigma       scaling parameter mutation distribution (0.1)
  ovote       'ignore' or 'account' ('account'); if set to 'ignore' group formation not dependent on x, not addressed in Port et al.
  bvote       one of 'ignore', 'kin', 'despotic', 'egalitarian' 'hierarchical' ('despotic'), if set to 'ignore', group formation not dependent on y (offspring control), if set to 'kin', acceptance of offspring to group only dependent on mother's y (mother's control), if set to 'egalitarian' acceptance of offspring dependent on average y of all group members (egalitarian group control), if set to 'hierarchical' acceptance of offspring dependent on female n's y, provided that female n-1 accepted offspring (hierarchical group control). Set to 'despotic' in Port et al.; other control systems not considered
  oplacement  offspring hierarchies placement 'back' or 'sort' ('sort'); if set to 'sort' offspring enters rank-hierarchy below mother's rank. Set to 'back' in Port et al.
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
