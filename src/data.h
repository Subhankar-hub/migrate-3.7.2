#ifndef DATA_INCLUDE
#define DATA_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 D A T A   R O U T I N E S 
 
 creates data structures,
 read data (Electrophoretic loci, sequences, microsats),
 prints data,
 destroys data.
 
 
 Theta(i)=4 N(i)mu
 M(ji) = m(ji)/mu
 4 N(i)m(ji) = Theta(i)M(ji)
                                                                                                               
*-----------------------------------------------------------------
  Bayesian and Maximum likelihood estimation of migration rates 
  using coalescent trees
 
  Peter Beerli
  Department of Scientific Computing
  Florida State University
  Tallahassee, FL 32306-4120
  beerli@fsu.edu
 
  
Copyright 1997-2002 Peter Beerli and Joseph Felsenstein
Copyright 2003-2008 Peter Beerli
Copyright 2009-2012 Peter Beerli and Michal Palczewski

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

$Id: data.h 2126 2013-01-09 03:33:51Z beerli $
 
*----------------------------------------------------------------
*/

#include "migration.h"

void create_data (data_fmt ** data);
void get_data (FILE * infile, data_fmt * data, option_fmt * options, world_fmt *world);
void print_data (world_fmt * world, option_fmt * options,
                     data_fmt * data);
void print_spectra(world_fmt * world, option_fmt * options,data_fmt * data);
void print_data_summary (FILE * file, world_fmt * world,
                             option_fmt * options, data_fmt * data);
long find_missing(data_fmt *data, long pop, long locus);
short findAllele (data_fmt * data, char s[], long locus);
void free_datapart (data_fmt * data, option_fmt * options, long locus);
void read_distance_fromfile (FILE * dfile, long tips, long nmlength,
                                 MYREAL **m);
void read_geofile (data_fmt * data, option_fmt * options, long numpop);

void init_data_structure1 (data_fmt ** data, option_fmt * options);
void init_data_structure2 (data_fmt ** data, option_fmt * options,
                               long pop);
void init_data_structure3 (data_fmt * data, option_fmt *options);
MYREAL create_alleles (data_fmt * data, option_fmt * options);
void set_numind (data_fmt * data);

long max_shuffled_individuals(option_fmt *options, data_fmt *data, long pop, long locus);

long number_genomes (int type);
void destroy_data(data_fmt * data);

#endif /*DATA_INCLUDE */
