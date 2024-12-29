/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 P R O F I L E    L I K E L I H O O D    R O U T I N E S 
 
 Peter Beerli 1997, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

$Id: profile.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/
extern boolean print_profile_likelihood_driver (long which, world_fmt * world,
            long *gmaxptr);
extern void print_profile_likelihood (long which, world_fmt * world, long *gmaxptr);
extern long print_profile_percentile (world_fmt * world);
extern void allocate_profile_percentiles (world_fmt * world);
extern void destroy_profile_percentiles (world_fmt * world);
extern void print_profile_title (world_fmt * world);
extern long warp (long ii);
extern long choose_profile_parameters (world_fmt * world, long *profilelist);

#define GRIDSIZE 9
#define GRIDMIDDLE 5
#define GRID    {0.01,0.05,0.10,0.50,0.99,0.95,0.90,0.50,1.0}
#define SHOWGRID {0.005,0.025,0.05,0.25,0.995,0.975,0.95,0.75,0.50}
#define GRID2   {0.005,0.025,0.05,0.25,0.5, 0.75,0.95,0.975,0.995}
#define INDEX {0,1,2,3,8,7,6,5,4}
#define DEVIATE {0.02,0.10,0.20, 0.5, 50.,10., 5., 2., 1.}
#define DEVIATE2 {0.02,0.10,0.20, 0.5, 1., 2., 5., 10., 50. }
//#define DEVIATE {0.002,0.010,0.020, 0.05, 5000.,1000., 500., 200., 1.}
//#define ABSOLUTE {1e-100,1e100}
