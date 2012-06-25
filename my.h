/*
  Calculate Mathematical Constants Using GMP (GNU Multiple Precision)

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

#ifndef _MY_H
#define _MY_H

#include "gmp.h"

#define  BITS_PER_DIGIT   3.32192809488736234787
#define  DOUBLE_PREC      53

void my_init(unsigned long prec);
void my_clear();

void my_sqrt_ui(mpf_t r, unsigned long x);
void my_invsqrt_ui(mpf_t r, unsigned long x);
void my_div(mpf_t r, mpf_t y, mpf_t x);
void my_sqrt(mpf_t r, mpf_t x);

void my_divexact(mpz_t r, mpz_t y, mpz_t x);
void my_div_set_threshold(unsigned long t);
unsigned long my_div_estimate_threshold();

void my_out_str(FILE *fp, unsigned long base, unsigned long digits, mpf_t f);

#endif

