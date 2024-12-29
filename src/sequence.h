#ifndef SEQUENCE_INCLUDE
#define SEQUENCE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 
 
 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
(c) 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
(c) 2003-2004 Peter Beerli, Tallahassee FL
 
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

$Id: sequence.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

extern void initratio (option_fmt * options);
extern void initfreqs (MYREAL *freqa, MYREAL *freqc, MYREAL *freqg,
                           MYREAL *freqt);
extern void initcatn (long *categs);
extern boolean initcategs (long categs, MYREAL *rate, MYREAL *probcat);
extern void initprobcat (long categs, MYREAL *probsum, MYREAL *probcat);
extern void constrain_rates(long categs, MYREAL *rate, MYREAL *probcat);
extern void initlambda (option_fmt * options);

extern void init_sequences (world_fmt * world, option_fmt * options,
			    data_fmt * data, long locus);
extern void init_sequences2 (world_fmt * world, seqmodel_fmt * seq,
			     long locus);
extern void init_tbl (world_fmt * world, long locus);
extern void print_weights (FILE * outfile, world_fmt * world,
			   option_fmt * options, long locus);
extern void print_tbl (FILE * outfile, world_fmt * world,
                           option_fmt * options, long locus);
extern void print_seqfreqs (FILE * outfile, world_fmt * world,
			    option_fmt * options);
extern MYREAL treelike_seq (world_fmt * world, long locus);
extern MYREAL treelike_snp (world_fmt * world, long locus);
extern void snp_invariants (contribarr invariants, world_fmt *world, long locus,  phenotype x1, MYREAL *scaler);
extern void make_sequences (world_fmt * world, option_fmt * options,
			    data_fmt * data, long locus);
extern void make_invarsites (world_fmt * world, data_fmt * data, long locus);
extern void make_invarsites_unlinked (world_fmt * world, data_fmt * data,
				      long locus);
extern void make_snp (world_fmt * world, option_fmt * options,
		      data_fmt * data, long locus);
extern MYREAL treelike_snp_unlinked (world_fmt * world, long locus);
extern void copy_seq (world_fmt * original, world_fmt * kopie);
extern void init_sequences_aliases (world_fmt *earth, world_fmt * world, option_fmt * options,
                                        data_fmt * data, long locus);
extern void find_rates_fromdata(data_fmt * data, option_fmt * options, world_fmt *world);
extern void find_rates_fromdata_alleles(data_fmt * data, option_fmt * options, world_fmt *world, MYREAL mean);
extern void free_seq(seqmodel_fmt **seq, long seqnum);

#endif /*SEQUENCE_INCLUDE */
