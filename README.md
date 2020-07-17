# npm
Natal philopatry model


## Building the simulation

If you are on Windows, you *can* skip the build-process and use the pre-compiled binary `npm.exe` in the `bin` folder. Otherwise, go ahead and build it with CMake:

```
:~/npm$ mkdir build && cd build
:~/npm/build$ cmake ..
:~/npm/build$ cmake --build . --config Release --target install
```

generating the Doxygen source code documentation (optional)
```
:~/npm$ doxygen .   # create documentation in ~/npm/doc/html
```

Tested on Linux (g++ > 8.0), MacOS (Xcode > 10) and Windows (Visual Studio 2019), this should have created the binary `:~npm/bin/npm`. If everything went well, you should be able to run:

```
:~/npm/bin$ ./npm --help
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
  phi         scramble competition parameter in F(n,R) (0.1)
  delta       contest competition parameter in F(n,R) (0.0)
  k           helping parameter in F(n,R) (10)
  Alleles     initial alleles A0,A1,A2,B0,B1,B2 as string ('5 0 0 5 0 0')
  Mask        alleles mask as string ('1 1 1 1 1 1')
  Sb          survival probability breeder (0.8)
  Sm          survival probability male (0.8)
  Sff         survival probability female floater (0.6)
  Smf         survival probability male floater (0.8)
  Smax        maximum survival (longevity) (0.95)
  sigma       shape parameter in Sx() (1.0)
  eps         patch search efficiency (0.005)
  t0          baseline takeover probability (0.05)
  tau         benefit for communal territory defense (1.0)
  mu          mutation probability (0.01)
  gamma       scaling parameter mutation distribution (0.1)
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
  npm -v mode=residency nmf=0 eps=0.0001 mu=0.001 gamma=0.01 ticks=1e6 log=10000 file=res2.R
```

## Result files

`npm ... file=res.R` creates an self-documented R-script, `res.R`, that can be directly sourced by R.

## Settings used in Port et al.

Please note that 'patches' are referred to 'home-ranges' in Port et al.

In the analyses presented in Port et al., the helping parameter k is always set to its default value (10), such that it has no effect on fecundity F(n,R). It is thus not described in equation (1). This parameter has been incorporated to allow for possible effects of alloparental care on breeder fecundity. Such effects are not considered in Port et al.

In addition to gene loci A0, A1, B0 and B1 mentioned in Port et al., our model consideres the effect of two additional gene loci: A2 and B2. A2 relates an established female's propensity to accept a juvenile, y(n), to the rank R of that juvenile's mother, and B2 relates a juvenile female's propensity to stay in the group, x(n), to her mother's rank. In the analyses presented in Port et al., the effect of A2 and B2 is masked by setting the Mask parameter to '1 1 0 1 1 0'. A2 und B2 are thus not described in equations (2a, b). Making the behavioral reaction norms dependent on rank R has only marginal effects on the results presented in Port et al.

In the analyses presented in Port et al., 'ovote' is always set to 'account' and 'bvote' is always set to 'despotic', representing the case where the decision over group membership depends on a juvenile female's propensity to stay in the group, x(n), and the most dominant female's propensity to accept the juvenile, y(n). Yet our model allows for other decision rules, too: By setting 'ovote' to 'ignore', the decision over group-membership depends solely on y(n). By setting 'bvote' to 'ignore', the decision over group-mebership depends solely on the juvenile's propensity to stay x(n) (offspring control over group-membership). By setting 'bvote' to 'kin' y(n) is the y(n) of the respective offspring's mother (rather than the y(n) of the most dominant female). By setting 'bvote' to 'egalitarian', y(n) is calculated as the average y(n) of all group members (group control). By setting 'bvote' to 'hierarchical', the decision over group-membership is made by all group member, but in a hierarchical way: y(n) is the y(n) of a female at rank R, given that the female at rank R-1 accepted the offspring. Different decision rules with respect to y(n) have only marginal effects of the results presented in Port et al.

In the analyses presented in Port et al., 'oplacemet' is always set to 'back', representing the case where new group members enter the group at the lowest rank position. Setting 'oplacement' to 'sort' means that new group members enter the group at the dominance rank directely below their mother. This model variant, too, has only marginal effects on the results presented in Port et al.

