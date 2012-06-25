/*
  Arithmetic with numbers in factorized form

  Copyright (C) 2006-2012 Hanhong Xue (macroxue at yahoo dot com)

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 59 Temple
  Place, Suite 330, Boston, MA 02111-1307 USA
*/
#ifndef _FAC_H
#define _FAC_H

#include "gmp.h"

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

typedef struct {
    unsigned max_facs;
    unsigned num_facs;
    unsigned *fac;
    unsigned *pow;
} fac_struct;

typedef fac_struct fac_t[1];

void fac_init_size(fac_t f, unsigned s);
void fac_init(fac_t f);
void fac_clear(fac_t f);

void fac_set_bp(fac_t f, unsigned base, unsigned pow);
void fac_mul_bp(fac_t f, unsigned base, unsigned pow);
void fac_mul(fac_t f, fac_t g);

void fac_gcd_compact(fac_t gcd, fac_t f, fac_t g);
void fac_to_mpz(mpz_t r, fac_t f);

void fac_show(fac_t f);

void fac_sieve_init(unsigned size);
void fac_sieve_clear();

#endif

