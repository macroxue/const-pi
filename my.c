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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "my.h"

mpf_t t1, t2;
mpf_t d1, d2;

void my_init(unsigned long prec)
{
    mpf_init2(t1, prec);
    mpf_init2(t2, prec);
    mpf_init2(d1, prec);
    mpf_init2(d2, prec);
}
void my_clear()
{
    mpf_clear(t1);
    mpf_clear(t2);
    mpf_clear(d1);
    mpf_clear(d2);
}

// r = sqrt(x)
void my_sqrt_ui(mpf_t r, unsigned long x)
{
    unsigned long prec, bits, prec0;

    prec0 = mpf_get_prec(r);

    if (prec0<=DOUBLE_PREC) {
        mpf_set_d(r, sqrt(x));
        return;
    }

    bits = 0;
    for (prec=prec0; prec>DOUBLE_PREC;) {
        int bit = prec&1;
        prec = (prec+bit)/2;
        bits = bits*2+bit;
    }

    mpf_set_prec_raw(t1, DOUBLE_PREC);
    mpf_set_d(t1, 1/sqrt(x));

    while (prec<prec0) {
        prec *=2;
        if (prec<prec0) {
            /* t1 = t1+t1*(1-x*t1*t1)/2; */
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, t1, t1);         // half x half -> full
            mpf_mul_ui(t2, t2, x);
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_div_2exp(t2, t2, 1);
            mpf_mul(t2, t2, t1);         // half x half -> half
            mpf_set_prec_raw(t1, prec);
            mpf_add(t1, t1, t2);
        } else {
            prec = prec0;
            /* t2=x*t1, t1 = t2+t1*(x-t2*t2)/2; */
            mpf_set_prec_raw(t2, prec/2);
            mpf_mul_ui(t2, t1, x);
            mpf_mul(r, t2, t2);          // half x half -> full
            mpf_ui_sub(r, x, r);
            mpf_mul(t1, t1, r);          // half x half -> half
            mpf_div_2exp(t1, t1, 1);
            mpf_add(r, t1, t2);
            break;
        }
        prec -= (bits&1);
        bits /=2;
    }
}

// r = 1/sqrt(x)
void my_invsqrt_ui(mpf_t r, unsigned long x)
{
    unsigned long prec, bits, prec0;

    prec0 = mpf_get_prec(r);

    if (prec0<=DOUBLE_PREC) {
        mpf_set_d(r, sqrt(x));
        return;
    }

    bits = 0;
    for (prec=prec0; prec>DOUBLE_PREC;) {
        int bit = prec&1;
        prec = (prec+bit)/2;
        bits = bits*2+bit;
    }

    mpf_set_prec_raw(t1, DOUBLE_PREC);
    mpf_set_d(t1, 1/sqrt(x));

    while (prec<prec0) {
        prec *=2;
        if (prec<prec0) {
            /* t1 = t1+t1*(1-x*t1*t1)/2; */
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, t1, t1);         // half x half -> full
            mpf_mul_ui(t2, t2, x);
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_div_2exp(t2, t2, 1);
            mpf_mul(t2, t2, t1);         // half x half -> half
            mpf_set_prec_raw(t1, prec);
            mpf_add(t1, t1, t2);
        } else {
            prec = prec0;
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, t1, t1);         // half x half -> full
            mpf_mul_ui(t2, t2, x);
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_div_2exp(t2, t2, 1);
            mpf_mul(t2, t2, t1);         // half x half -> half
            mpf_add(r, t1, t2);
            break;
        }
        prec -= (bits&1);
        bits /=2;
    }
}

// r = y/x   WARNING: r cannot be the same as y.
void my_div(mpf_t r, mpf_t y, mpf_t x)
{
    unsigned long prec, bits, prec0;
    assert(r != y);

    prec0 = mpf_get_prec(r);

    if (prec0<=DOUBLE_PREC) {
        mpf_set_d(r, mpf_get_d(y)/mpf_get_d(x));
        return;
    }

    bits = 0;
    for (prec=prec0; prec>DOUBLE_PREC;) {
        int bit = prec&1;
        prec = (prec+bit)/2;
        bits = bits*2+bit;
    }

    mpf_set_prec_raw(t1, DOUBLE_PREC);
    mpf_set_d(t1, 1/mpf_get_d(x));

    while (prec<prec0) {
        prec *=2;
        if (prec<prec0) {
            /* t1 = t1+t1*(1-x*t1); */
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, x, t1);          // full x half -> full
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_mul(t2, t2, t1);         // half x half -> half
            mpf_set_prec_raw(t1, prec);
            mpf_add(t1, t1, t2);
        } else {
            prec = prec0;
            /* t2=y*t1, t1 = t2+t1*(y-x*t2); */
            mpf_set_prec_raw(t2, prec/2);
            mpf_mul(t2, t1, y);          // half x half -> half
            mpf_mul(r, x, t2);           // full x half -> full
            mpf_sub(r, y, r);
            mpf_mul(t1, t1, r);          // half x half -> half
            mpf_add(r, t1, t2);
            break;
        }
        prec -= (bits&1);
        bits /=2;
    }
}

// r = sqrt(x)
void my_sqrt(mpf_t r, mpf_t x)
{
    unsigned prec, bits, prec0;

    prec0 = mpf_get_prec(r);

    if (prec0 <= DOUBLE_PREC) {
        mpf_set_d(r, sqrt(mpf_get_d(x)));
        return;
    }

    bits = 0;
    for (prec = prec0; prec > DOUBLE_PREC;) {
        int bit = prec & 1;
        prec = (prec + bit) / 2;
        bits = bits * 2 + bit;
    }

    mpf_set_prec_raw(t1, DOUBLE_PREC);
    mpf_set_d(t1, 1 / sqrt(mpf_get_d(x)));

    while (prec < prec0) {
        prec *= 2;
        /*printf("prec=%d, prec0=%d\n", prec, prec0); */
        if (prec < prec0) {
            /* t1 = t1+t1*(1-x*t1*t1)/2; */
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, t1, t1);
            mpf_set_prec_raw(x, prec/2);
            mpf_mul(t2, t2, x);
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_div_2exp(t2, t2, 1);
            mpf_mul(t2, t2, t1);
            mpf_set_prec_raw(t1, prec);
            mpf_add(t1, t1, t2);
        } else {
            prec = prec0;
            /* t2=x*t1, t1 = t2+t1*(x-t2*t2)/2; */
            mpf_set_prec_raw(t2, prec/2);
            mpf_set_prec_raw(x, prec/2);
            mpf_mul(t2, t1, x);
            mpf_mul(r, t2, t2);
            mpf_set_prec_raw(x, prec);
            mpf_sub(r, x, r);
            mpf_mul(t1, t1, r);
            mpf_div_2exp(t1, t1, 1);
            mpf_add(r, t1, t2);
            break;
        }
        prec -= (bits & 1);
        bits /= 2;
    }
}

// Exact Division
unsigned long div_threshold = 1404000;
void my_divexact(mpz_t r, mpz_t y, mpz_t x)
{
    unsigned long prec = mpz_sizeinbase(y, 2);
    unsigned long prec2 = mpz_sizeinbase(x, 2);
    if (prec >= div_threshold) {
        mpf_set_prec_raw(d1, prec);
        mpf_set_prec_raw(d2, prec);
        mpf_set_z(d1, y);
        mpf_set_z(d2, x);
        mpf_div_2exp(d1, d1, prec);
        mpf_div_2exp(d2, d2, prec2);
        my_div(d2, d1, d2);
        mpf_mul_2exp(d2, d2, prec-prec2);
        mpf_set_d(d1, 0.5);
        mpf_add(d2, d2, d1);
        mpz_set_f(y, d2);
    } else {
        mpz_divexact(y, y, x);
        //mpz_tdiv_q(y, y, x);
    }
}

unsigned long my_div_estimate_threshold()
{
    const double log10 = 3.321928;
    clock_t begin, end, diff1, diff2;
    mpz_t y, x;
    unsigned long digits = 100, min_digits, max_digits, prec, stage = 1;
    div_threshold = 0;
    for (;;) {
        if (stage == 1)
            digits *= 2;
        else
            digits = (min_digits + max_digits)/2;
        prec = (unsigned long) (digits * log10) + 32;
        my_init(prec);
        mpz_init(y);
        mpz_init(x);
        mpz_ui_pow_ui(y, 10, digits);
        mpz_ui_pow_ui(x, 10, digits/10);
        begin = clock();
        my_divexact(y, y, x);
        end = clock();
        diff1 = end - begin;

        mpz_ui_pow_ui(y, 10, digits);
        begin = clock();
        mpz_divexact(y, y, x);
        end = clock();
        diff2 = end - begin;
        my_clear();
        mpz_clear(y);
        mpz_clear(x);

        printf("digits %ld, my %ld, mpz %ld\n", digits, diff1, diff2);
        if (stage == 1) {
            if (diff1 >= 10 && diff1 < diff2) {
                stage = 2;
                min_digits = digits/2;
                max_digits = digits;
            }
        } else {
            if (diff1 < diff2)
                max_digits = digits;
            else
                min_digits = digits;
            if (max_digits - min_digits < 100)
                break;
        }
    }
    return (unsigned long) (digits * log10);
}

void my_div_set_threshold(unsigned long t)
{
    div_threshold = t;
}

// Convert Hexadecimal to Decimal
#define UNIT_SIZE    5
#define UNIT_MOD     (int)1e5
#define LINE_SIZE    (UNIT_SIZE * 10)
#define NUM_BLOCKS   4
char   out_buf[(LINE_SIZE+30)*10], *out_ptr = out_buf;
mpf_t  trunk;

void utoa(unsigned long i, unsigned long d)
{
    static unsigned long pow10[] = { 0, 10, 100, 1000, 10000, 100000,
                                     1000000, 10000000, 100000000, 1000000000
                                   };
    if (d == 0) {
        for (d = 0; d < sizeof(pow10)/sizeof(pow10[0]); d++)
            if (i < pow10[d]) break;
    }
    while (--d) {
        *out_ptr = i / pow10[d];
        i = i - *out_ptr * pow10[d];
        *out_ptr++ += '0';
    }
    *out_ptr++ = i + '0';
}

void flush_out(FILE *fp)
{
    *out_ptr++ = '\n';
    *out_ptr = '\0';
    fputs(out_buf, fp);
    out_ptr = out_buf;
}

void my_out_str_raw(FILE *fp, unsigned long digits, mpf_t f, unsigned long offset)
{
    unsigned long d;

    if (digits <= LINE_SIZE*NUM_BLOCKS) {
        unsigned long cursor = offset % LINE_SIZE;
        for (d = 0; d < digits; ) {
            mpf_set_prec_raw(f, (int)((digits-d)*BITS_PER_DIGIT+1));
            mpf_mul_ui(f, f, UNIT_MOD);
            unsigned long i = mpf_get_ui(f);
            mpf_sub_ui(f, f, i);

            utoa(i, UNIT_SIZE);
            *out_ptr++ = ' ';
            d += UNIT_SIZE;
            cursor += UNIT_SIZE;
            if (cursor == LINE_SIZE) {
                cursor = 0;
                *out_ptr++ = ':';
                *out_ptr++ = ' ';
                utoa(offset + d, 0);
                *out_ptr++ = '\n';
                if ((offset + d) % (LINE_SIZE*10) == 0)
                    flush_out(fp);
            }
        }
    } else {
        mpf_t block, mod;
        unsigned long num_units = (digits + UNIT_SIZE-1)/UNIT_SIZE;
        unsigned long block_size =  (num_units + NUM_BLOCKS-1)/NUM_BLOCKS*UNIT_SIZE;
        mpf_set_default_prec((int)(block_size*BITS_PER_DIGIT+1));
        mpf_init(block);
        mpf_init_set_ui(mod, 10);
        mpf_pow_ui(mod, mod, block_size);

        for (d = 0; d < digits; d += block_size) {
            unsigned long size = block_size < digits - d ? block_size : digits - d;
            mpf_set_prec_raw(block, (int)(size*BITS_PER_DIGIT+1));
            mpf_set(block, f);
            my_out_str_raw(fp, size, block, offset+d);
            if (block_size < digits - d) {
                mpf_set_prec_raw(f, (int)((digits-d)*BITS_PER_DIGIT+1));
                mpf_mul(f, f, mod);
                mpf_floor(trunk, f);
                mpf_sub(f, f, trunk);
            }
        }
        mpf_clear(block);
        mpf_clear(mod);
    }
}

void my_out_str(FILE *fp, unsigned long base, unsigned long digits, mpf_t f)
{
    unsigned long num_units = (digits + UNIT_SIZE-1)/UNIT_SIZE;
    unsigned long block_size =  (num_units + NUM_BLOCKS-1)/NUM_BLOCKS*UNIT_SIZE;
    mpf_init2(trunk, (int)(block_size*BITS_PER_DIGIT+1));

    unsigned long i = mpf_get_ui(f);
    fprintf(fp, "%lu.\n", i);
    mpf_sub_ui(f, f, i);
    my_out_str_raw(fp, digits, f, 0);
    if (out_ptr != out_buf)
        flush_out(fp);

    mpf_clear(trunk);
}
