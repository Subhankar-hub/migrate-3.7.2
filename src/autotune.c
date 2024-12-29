/*
 Bayesian inference or maximum likelihood inference of structured population models
 
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
 
 
 */

#include "migration.h"
#include "tools.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#include "bayes.h"
#include "mcmc.h"
#include "reporter.h"

extern void print_menu_equilib (world_fmt * world);
extern int myID;

void autotune_proposal(world_fmt *world, long which);
void present_burnin_info(world_fmt *world, MYREAL ess, MYREAL acceptance, MYREAL var, MYREAL oldvar, long step);
void burnin_ml(world_fmt * world);
void burnin_progress(world_fmt *world, MYREAL percent);
void burnin_bayes(world_fmt * world);
MYREAL mean_acceptance_rate(world_fmt * world);
boolean auto_stop_burnin(world_fmt *world, long step, long stop, MYREAL * var, MYREAL *autocorrelation, MYREAL * effective_sample, MYREAL *acceptances);
void burnin_chain (world_fmt * world);


void autotune_proposal(world_fmt *world, long which)
{
    MYREAL ratio;
    worldoption_fmt *wopt = world->options;
    bayes_fmt *bayes = world->bayes;
    MYREAL *delta = bayes->delta;
    long space = MAX(10,world->numpop2+1);
    const MYREAL ma = bayes->maxparam[which];
    const MYREAL mi = bayes->minparam[which];
    const MYREAL mindelta = (ma-mi)/1000.;
    MYREAL olddelta;
    if(world->in_burnin && world->cold)
    {
        // formula: delta = (pR - 1) * (min-max)
        // for pR=0.44 and min=0 and max=0.1 --> delta= -0.66 * -0.1 = 0.066
        // pR = delta/(min-max) + 1
        // Binom(n,k) pR^k (1-pR)^(n-k)
        // n=1:k=1: delta/(min-max)+1
        //      delta/(-10*delta)+1 ==> 1-1/10 = 0.9
        // n=1:k=0: (1-delta/(min-max)-1)) ==> -1/10 = -0.1
        if(wopt->has_autotune)
        {
	  if(bayes->trials[which] % space == 0)
	    {
	      olddelta = delta[which];
	      ratio = bayes->accept[which]/(1.0+bayes->trials[which]);
	      if(wopt->autotune < ratio)
		{
		  delta[which] *= 1.0101; /*0.99^(r/(r-1))*/
		}
	      else
		{
		  delta[which] *= 0.990;
		}
	      if(delta[which] > ma)
                delta[which] = ma;
	      if(delta[which] < mindelta)
		delta[which] = olddelta;
	    }
	}
    }
}

void present_burnin_info(world_fmt *world, MYREAL ess, MYREAL acceptance, MYREAL var, MYREAL oldvar, long step)
{
  long z=0;
#ifdef MPI
    char p[LINESIZE];
    char p1[LINESIZE];
    long bufsize;
    if(myID!=MASTER)
    {
      bufsize = sprintf(p,"%li %f %f %f %f %li\n",world->locus, ess, acceptance, var,oldvar,step);
        sprintf(p1,"B%li",bufsize);
        MYMPISEND (p1, SMALLBUFSIZE, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+BURNTAG, comm_world);
        MYMPISEND (p, bufsize, MPI_CHAR, (MYINT) MASTER, (MYINT) myID+BURNTAG, comm_world);
    }
#endif
    z = world->burnin_z;
    world->burnin_stops[z].locus     =  world->locus;
    world->burnin_stops[z].replicate =  world->rep;
    world->burnin_stops[z].stopstep  = step;
    world->burnin_stops[z].ess       = ess;
    world->burnin_stops[z].accept  = acceptance;
    world->burnin_stops[z].variance  = var;
    world->burnin_stops[z].oldvariance = oldvar;
    world->burnin_stops[z].worker    = myID;
    world->burnin_z++;
}

MYREAL mean_acceptance_rate(world_fmt * world)
{
  long which;
  long i;
  bayes_fmt * bayes = world->bayes;
  MYREAL sum = 0.0;
  long count = 0;
  if (!world->options->bayes_infer)
    return -1.0;

  for (i = 0; i <  world->numpop2 + bayes->mu; i++)
    {
      if(shortcut(i,bayes,&which))
	continue;
      else
	{
	  sum += world->bayes->accept[which]/(1.0+world->bayes->trials[which]);
	  count += 1;
	}
    }
  return sum/count;
}

boolean auto_stop_burnin(world_fmt *world, long step, long stop, MYREAL * var, MYREAL *autocorrelation, MYREAL * effective_sample, MYREAL *acceptances)
{
  (void) acceptances;
    char autostop = world->options->burnin_autostop;
    const long delta = ((stop >= 10000) ? 1000 : (stop / 10));
    const long nn = world->numpop2 + world->bayes->mu + 1 ;
    MYREAL oldvar= *var;
    MYREAL acceptance;
    MYREAL ess;
    boolean done = FALSE;
    boolean acceptanceOK = FALSE;
    boolean essOK = FALSE;
    boolean vardiffOK = FALSE;
    single_chain_var (world, step, var, autocorrelation, effective_sample);
    if (step > delta)
      vardiffOK = (fabs(*var/oldvar - 1.0) < world->varheat);
    else
      vardiffOK = FALSE;
    essOK = max_ess(effective_sample,nn,world->essminimum, &ess);
    acceptance = mean_acceptance_rate(world);
    acceptanceOK = acceptance > world->options->autotune - 0.05 && acceptance < world->options->autotune + 0.05;  
    switch(autostop)
    {
    case 'a':
      done = vardiffOK;
      break;
    case 't':
      done = acceptanceOK;
      break;
    case 'e':
      done = essOK;
      break;
    case ' ':
    default:
      break;
    }
    if(done || step==stop)
    {
      if(world->cold)
	present_burnin_info(world,ess, acceptance, *var,oldvar,step);
      return done;
    }
    return FALSE;
}

void burnin_ml(world_fmt * world)
{
  boolean done=FALSE;
  MYREAL *autocorrelation;
  MYREAL *effective_sample;
  MYREAL *acceptances;
  MYREAL var=1000000.0;
  long step;
  const long stop  = world->options->burn_in * world->increment;
  long nn = world->numpop2 + world->bayes->mu;
  long delta = ((stop > 10000) ? 1000 : (stop / 10));
  if(stop==0)
    return;
  if(delta==0)
    delta=1;

  autocorrelation  = (MYREAL *) mycalloc((3L*nn),sizeof(MYREAL));
  effective_sample = autocorrelation + nn;
  acceptances      = effective_sample + nn;

  single_chain_var (NULL, 0, &var, NULL, NULL);
  for(step=0;step<stop+1; step++)
    {
      //includes autotuning
      tree_update (world, 0);
      if((step % delta) == 0 && world->cold)
        {
	  done = auto_stop_burnin(world, step, stop, &var, autocorrelation, effective_sample, acceptances);  
        }
      if(done)
	break;
    }
  myfree(autocorrelation);
}

void burnin_progress(world_fmt *world, MYREAL percent)
{
  (void) world;
    char nowstr[STRSIZE]; // will hold time of day
    get_time (nowstr, "%H:%M:%S");
#ifdef MPI
    FPRINTF(stdout, "[%3i] %s   Burn-in %3.1f%% complete\n",myID, nowstr, percent); 
#else
    FPRINTF(stdout, "%s   Burn-in %3.1f%% complete\n",nowstr, percent); 
#endif
}


void burnin_bayes(world_fmt * world)
{
#ifndef MPI
  const boolean progress = world->options->progress && world->cold;
#endif
  const long stop = world->options->burn_in * world->increment;
  const long nn = world->numpop2 + world->bayes->mu;
  long delta = ((stop > 100000) ? (stop / 10) : (stop/3));
  MYREAL * autocorrelation;
  MYREAL * effective_sample;
  MYREAL * acceptances;
  MYREAL var= (MYREAL) ((MYREAL) HUGE);
  long step;
  boolean done=FALSE;

  if(stop==0)
    return;
  if(delta==0)
    delta=1;

  autocorrelation  = (MYREAL *) mycalloc((3L + 3L*nn),sizeof(MYREAL));
  effective_sample = autocorrelation + nn + 1;
  acceptances      = effective_sample + nn + 1;

  single_chain_var (NULL, 0, &var, NULL, NULL);
  for(step=0;step<stop; step++) // Cesky Krumlov 2013
    {
#ifndef MPI
      // too much writing for MPI therefore we do not report progress on burnin
      if (progress && (step+1) % delta == 0)
	burnin_progress(world,100*step/stop);
#endif
      //includes autotuning
        world->bayesaccept = bayes_update (world);
        if(world->bayesaccept == -1) //either do tree updates or parameter updates
        {
            tree_update (world, 0);
        }
        if((step % delta) == 0)
        {
	  done = auto_stop_burnin(world, step, stop, &var, autocorrelation, effective_sample, acceptances);
        }
	if(done)
	  break;
        world->bayes->count = 0;
    }
  myfree(autocorrelation);
}

void burnin_chain (world_fmt * world)
{
  char nowstr[STRSIZE];
  char autostop = world->options->burnin_autostop;
  size_t nng = (size_t) (world->numpop2 + world->bayes->mu + 1); //to set zero including the genealogy
  const boolean treeprint = (boolean) world->options->treeprint;
    long z=0;
    world->burnin_z=0;        
    world->options->treeprint = NONE;
    
    world->in_burnin = TRUE;
    if (world->cold)
    {
        print_menu_equilib (world);
    }
    if(world->options->bayes_infer)
    {
        burnin_bayes(world);
	memset(world->bayes->accept,0,sizeof(long)*nng);//Cesky Krumlov 2013
	memset(world->bayes->trials,0,sizeof(long)*nng);//Cesky Krumlov 2013
    }
    else
    {
      burnin_ml(world);
    }
    if (world->cold)
      {
	z = world->burnin_z-1;
	switch(autostop)
	  {
	  case 'a':
	    FPRINTF(stdout,"[%3i]            stopped at step %li with var-ratio=%.2f/%.2f=%.3f (min(ESS)=%.2f, avg(acceptance)=%.2f)\n",
		    myID,world->burnin_stops[z].stopstep,world->burnin_stops[z].variance,world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].ess, world->burnin_stops[z].accept);
	    break;
	  case 't':
	    FPRINTF(stdout,"[%3i]            stopped at step %li with avg(acceptance)=%.2f (var-ratio=%.2f/%.2f=%.3f, min(ESS)=%.2f)\n",
		    myID,world->burnin_stops[z].stopstep,
		    world->burnin_stops[z].accept,
		    world->burnin_stops[z].variance,world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].ess);
	    break;
	  case 'e':
	    FPRINTF(stdout,"[%3i]            stopped at step %li with min(ess)=%.2f (var-ratio=%.2f/%.2f=%.3f, avg(acceptance)=%.2f)\n",
		    myID,world->burnin_stops[z].stopstep,
		    world->burnin_stops[z].ess,
		    world->burnin_stops[z].variance,world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].variance/world->burnin_stops[z].oldvariance,
		    world->burnin_stops[z].accept);
	    break;
	  default:
	    break;
	  }
	get_time (nowstr, "%H:%M:%S");
#ifdef MPI
	FPRINTF(stdout,"[%3i] %8.8s   Sampling ...\n",myID, nowstr);
#else
	FPRINTF(stdout,"%8.8s   Sampling ...\n", nowstr);
#endif

      }
    world->options->treeprint = treeprint;
    world->in_burnin = FALSE;
}
