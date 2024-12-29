#ifndef _GAMMALIKE_H_
#define _GAMMALIKE_H_
/* gamma deviated mutation rate among loci


 code for the calculation of the gamma likelihood
 and for its derivatives
 
(c) Peter Beerli 2013 Tallahassee
 
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
*/

#include "migration.h"

MYREAL gamma_loci_like (helper_fmt * helper, MYREAL *oparam,
                                   MYREAL *olparam, MYREAL denom);
MYREAL gamma_locus_like (nr_fmt * nr, MYREAL *oparam, MYREAL *olparam,
                                    MYREAL denom, long locus);
void gamma_loci_derivative (helper_fmt * helper);
void gamma_locus_derivative (helper_fmt * helper, long locus);
void gamma_loci_difference (helper_fmt * helper);
#endif
