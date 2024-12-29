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
/*! \file priors.c 

priors.c supplies prior distributions
*/
#include "migration.h"
#include "random.h"
#include "definitions.h"
#include "tools.h" 

#ifdef HAVE_LGAMMA
#define LGAMMA lgamma
#else
#define LGAMMA mylgamma
#endif

///
/// Normal distribution
/// this generates only one random number, but could be easily changed to supply two numbers
MYREAL normal_rand(MYREAL mean, MYREAL std)
{
  MYREAL u = RANDUM();
  MYREAL v = RANDUM();
  return sqrt(-2 * log(u)) * cos(TWOPI*v);
  //return sqrt(-2 * log(u)) * sin(TWOPI*v);
}

///
// Gamma distribution
//
// gamma deviated random number 
// from Marsaglia G., Tsang, W. W. (2000) A simple method for generating Gamma variables. 
// ACM Transactions on Mathematical Software Vol 26. No. 3: 363-372
// requirements are a 
//   - normal random number generator
//   - uniform random number generator
//
// [see Python program to test implementation] 
//
// returns a random number for shape parameter alpha>=1 and scale parameter beta=1
MYREAL __gamma_rand(MYREAL alpha)
{
  const MYREAL d = alpha - 1./3.;
  const MYREAL c = 1./sqrt(9. * d);
  MYREAL v,x,xx, u;
  while(1)
    {
      x = normal_rand(0.,1.);
      v = 1. + c * x;
      while(v<=0.0)
	{
	  x = normal_rand(0.,1.);
	  v = 1. + c * x;
	}
      v  = v * v * v;
      u  = RANDUM();
      xx = x * x;
      if (u < 1.0-0.0331*xx*xx)
	return d*v;
      if (log(u) < 0.5*xx + d * (1.0 - v + log(v)))
	return d*v;
    }
}

/// returns a random number from a gamma distribution with shape parameter alpha and scale parameter beta
MYREAL gamma_rand(MYREAL alpha, MYREAL beta)
{
  MYREAL aa;
  if (alpha<1.0)
    {
      aa = 1.0 + alpha;
      return __gamma_rand(aa)*pow(RANDUM(),(1.0/alpha)) * beta;
    }
  else
    {
      return __gamma_rand(alpha) * beta;
    }
}

/// returns a random number from a truncated gamma distribution
MYREAL trunc_gamma_rand(MYREAL alpha, MYREAL beta, MYREAL lower, MYREAL upper)
{
  MYREAL x;
  while(1)
    {
      x = gamma_rand(alpha,beta);
      if (lower < x)
	{
	  if (x < upper)
	    return x;
	}
    }
}

MYREAL logpdf_gamma(MYREAL a, MYREAL b, MYREAL x)
{
  if (a>0.0 && b> 0.0)
    return -(x/b) - a*log(b) + (-1 + a)*log(x) - LGAMMA(a);
  else 
    return -((MYREAL) HUGE);
}

MYREAL cdf_gamma(MYREAL a, MYREAL b, MYREAL x)
{
  //Gamma[a, 0, xmax/b]/Gamma[a]
  return incompletegamma(x/b,a);
}

MYREAL logpdf_truncgamma(MYREAL a, MYREAL b, MYREAL  xmin, MYREAL xmax, MYREAL x)
{
  if((x > xmin) && (x < xmax))
    return logpdf_gamma(a,b,x)-log(cdf_gamma(a,b,xmax)-cdf_gamma(a,b,xmin));
  else 
    return -((MYREAL) HUGE);
}


MYREAL mean_truncgamma(MYREAL alpha, MYREAL beta, MYREAL lower, MYREAL upper)
{
  MYREAL m;
  MYREAL nom;
  MYREAL denom ;
  //1>  beta*(Gamma(1+alpha,lower/beta,upper/beta)/Gamma(alpha,lower/beta,upper/beta)
  //2> (beta*(Gamma(1 + alpha,upper/beta) - Gamma(1 + alpha,lower/beta)))/
  //   (Gamma(alpha,upper/beta) - Gamma(alpha,lower/beta));
  // 1 and 2 should be the same
  // 
  nom = beta * alpha *(incompletegamma(upper/beta,1+alpha) - incompletegamma(lower/beta,1+alpha));
  denom =  incompletegamma(upper/beta,alpha) - incompletegamma(lower/beta,alpha);
  m = nom / denom;
  return m;
}


MYREAL find_beta_truncgamma(MYREAL mean, MYREAL alpha, MYREAL lower, MYREAL upper)
{
  /* from wikipedia:
     INPUT: Function f, endpoint values a, b, tolerance TOL, maximum iterations NMAX
     CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0
     OUTPUT: value which differs from a root of f(x)=0 by less than TOL
     N ← 1
     While N ≤ NMAX { limit iterations to prevent infinite loop
     c ← (a + b)/2 new midpoint
     If (f(c) = 0 or (b – a)/2 < TOL then { solution found
     Output(c)
     Stop
     }
     N ← N + 1 increment step counter
     If sign(f(c)) = sign(f(a)) then a ← c else b ← c new interval
     }
     Output("Method failed.") max number of steps exceeded  
  */
  MYREAL tolerance = EPSILON;
  MYREAL nmax = 100;
  long n=1;
  MYREAL beta; /*=c in readme*/
  MYREAL fbeta;
  MYREAL fa;
  MYREAL a = lower;
  MYREAL b = upper;
  
  if (a==0.0)
    a += SMALLEPSILON;
  while (n<= nmax)
    {
      beta = (a+b)/2.0;
      fbeta = mean_truncgamma(alpha,beta,lower,upper) - mean;    
      if(fbeta  == 0.0 || (b-a)/2. < tolerance)
	return beta;
      n += 1;
      fa = mean_truncgamma(alpha,a,lower,upper) - mean;    
      if(((fbeta<0.0) && (fa<0.0)) || ((fbeta>0.0) && (fa>0.0)))
	a = beta;
      else
	b = beta;
    }
  return beta;
}
	

///
/// Gamma prior retrieve a new value from a truncated gamma between lower and upper
/// currently does not use the old parameter, but this 
MYREAL
propose_gamma_newparam (MYREAL param, long which, bayes_fmt * bayes, MYREAL *r)
{
  //const MYREAL delta = bayes->delta[which];
  const MYREAL lower = bayes->minparam[which];
  const MYREAL upper = bayes->maxparam[which];
  //const MYREAL mean = bayes->meanparam[which];
  MYREAL alpha = bayes->alphaparam[which];
  MYREAL beta = bayes->betaparam[which];
  return trunc_gamma_rand(alpha, beta, lower, upper);
}



///
/// Hastings ratio calculator for gamma distribution
/// P(new -> old)    P(old)
/// ------------- = -------
/// P(old -> new)    P(new)
/// cancels with log_prior_gamma -> 0.0
MYREAL hastings_ratio_gamma(MYREAL newparam, MYREAL oldparam, MYREAL delta, MYREAL r, bayes_fmt * bayes, long whichparam)
{
  if((newparam > bayes->maxparam[whichparam]) || (newparam < bayes->minparam[whichparam]))
    return -((MYREAL) HUGE);
  else
    return 0.;
}


///
/// Log Prior gamma distribution ratios between old and new parameter:
/// cancels with hastings ratio
MYREAL log_prior_ratio_gamma(MYREAL newparam, 
			   MYREAL oldparam, 
			   bayes_fmt * bayes, 
			   long which)
{
  if((newparam > bayes->maxparam[which]) || (newparam < bayes->minparam[which]))
    return -((MYREAL) HUGE);
  else
    return 0.;
}

///
/// Gamma prior distribution for theta or migration rate used in heating acceptance
/// uses pdf_truncgamma(a,b,min,max,x)
MYREAL log_prior_gamma(world_fmt *world, long numparam)
{
  //  long frompop;
  //long topop;
  MYREAL p0;
  long numpop = world->numpop;
  long start = ((numparam <= numpop || numpop==1) ? 0 : numpop);
  long stop = ((start == numpop) ? world->numpop2 : numpop);
  long i;
  MYREAL * param0 = world->param0;
  bayes_fmt * bayes = world->bayes;
  MYREAL a;
  MYREAL b;
  MYREAL val=0.0;
  for(i = start; i < stop; i++)
    {
      if(!strchr("0c", world->options->custm2[i]))
	{
	  //@@if(i>=numpop)//@@ && !world->options->usem)
	  //@@  {
#ifndef PRIORTEST 
	      //@@     m2mm(i,numpop,&frompop,&topop);
#endif
	      //@@p0 = param0[i] * param0[topop];
	  //@@ }
	  //@@else
	  //@@{
	      p0 = param0[i];
	      //@@  }

	  if((p0 > bayes->maxparam[i]) || (p0 < bayes->minparam[i]))
	    return -((MYREAL) HUGE);
	  else
	    {
	      a = bayes->alphaparam[i];
	      b = bayes->betaparam[i];
	      val += logpdf_truncgamma(a,b,bayes->minparam[i],bayes->maxparam[i],p0);
	      //old!! val += -p0 * ib + a * log(ib) + (a - 1.) * log(p0) - LGAMMA(a) ;
	    }
	}
    }
  return val;
}

///
/// 
MYREAL log_prior_gamma1(world_fmt *world, long numparam, MYREAL val)
{
  bayes_fmt * bayes = world->bayes;
  MYREAL retval;
  MYREAL a = bayes->alphaparam[numparam];
  MYREAL b = bayes->betaparam[numparam]; 
  retval =  logpdf_truncgamma(a,b,bayes->minparam[numparam],bayes->maxparam[numparam],val);
  return retval;
}

// test section
// compile using this:
// gcc -o priortest -g priors.c random.c sighandler.c tools.c -DPRIORTEST -DMERSENNE_TWISTER -DNOJPEGLIB -DMEXP=19937
// then call by  
// priortest alpha beta
// it will print 10000 numbers
// the mean and standard deviation
//
#ifdef PRIORTEST
#include <stdio.h>
#include "sighandler.h"

char * generator;
int myID;
long *seed;

 int main(long argc, char **argv)
 {
   long i;
   MYREAL xx;
   MYREAL pxx;
   MYREAL a;
   MYREAL b;
   MYREAL mean=0.0;
   MYREAL var=0.0;
   
   world_fmt *world;
   world = calloc(1,sizeof(world_fmt));
   world->bayes = calloc(1,sizeof(bayes_fmt));
   world->bayes->maxparam = calloc(1,sizeof(MYREAL));
   world->bayes->minparam = calloc(1,sizeof(MYREAL));
   world->bayes->alphaparam = calloc(1,sizeof(MYREAL));
   world->bayes->betaparam = calloc(1,sizeof(MYREAL));
   world->bayes->meanparam = calloc(1,sizeof(MYREAL));
   world->numpop=1;

   generator = (char *) mycalloc (1,sizeof(char) * 80);

   init_gen_rand(123789);

   a = atof(argv[1]);
   b = atof(argv[2]);
   world->bayes->alphaparam[0] = a;
   world->bayes->meanparam[0] = a*b; 
   world->bayes->minparam[0] = 1.0; 
   world->bayes->maxparam[0] = 50.0; 
   world->bayes->betaparam[0] = find_beta_truncgamma(a*b, a, world->bayes->minparam[0],world->bayes->maxparam[0]);
   printf("Truncated gamma distribution with alpha=%f, beta=%f, lower=%f, upper=%f\n",
	  a,world->bayes->betaparam[0],world->bayes->minparam[0],world->bayes->maxparam[0]);
   b = world->bayes->betaparam[0];
   for (i=0;i<10000;i++)
     {
       xx=trunc_gamma_rand(a,b,world->bayes->minparam[0], world->bayes->maxparam[0]);
       pxx = log_prior_gamma1(world, 0, xx);
       mean += (xx - mean)/(i+1);
       var  += (xx - mean) * (xx-mean);
       printf("%f %f\n",xx, pxx);
     }
   printf("results of random gamma truncated 0 and 500, using a=%f, b=%f\n",a,b);
   printf("Mean=%f Standard deviation = %f\n",mean, sqrt(var/10000));
   printf("Expected=%f expected standard deviation = %f\n",a*b, sqrt(a)*b);
   MYREAL v = (world->bayes->maxparam[0]-world->bayes->minparam[0]) * 0.2;
   MYREAL fx = cdf_gamma(a,b,0.5);
   MYREAL fx2 = logpdf_gamma(a,b,v);
   MYREAL fx3 = logpdf_truncgamma(a,b,world->bayes->minparam[0],world->bayes->maxparam[0],v);
   printf("v=%f CDF(%f,%f,%f)=%f PDF(%f,%f,%f)=%f TPDF(%f,%f,%f,%f,%f)=%f\n",v,a,b,0.5,fx,a,b,v,fx2,a,b,
	  world->bayes->minparam[0],world->bayes->maxparam[0],v,fx3);
 }
#endif
