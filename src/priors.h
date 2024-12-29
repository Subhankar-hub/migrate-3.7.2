#ifndef PRIORS_H
#define PRIORS_H
/*------------------------------------------------------
Bayesian inference 
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo algorithm
-------------------------------------------------------
    Prior distributions   R O U T I N E S

send questions concerning this software to:
Peter Beerli
beerli@fsu.edu

Copyright 2013 Peter Beerli 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 $Id$
    */
/*! file priors.c 

priors.c supplies prior distributions
*/

extern MYREAL find_beta_truncgamma(MYREAL mean, MYREAL alpha, MYREAL lower, MYREAL upper);
extern MYREAL propose_gamma_newparam (MYREAL param, long which, bayes_fmt * bayes, MYREAL *r);
extern MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
extern MYREAL log_prior_ratio_gamma(MYREAL newparam, MYREAL oldparam,  bayes_fmt * bayes,  long which);
extern MYREAL log_prior_gamma(world_fmt *world, long numparam);
extern MYREAL log_prior_gamma1(world_fmt *world, long numparam, MYREAL val);
extern MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam);
#endif
