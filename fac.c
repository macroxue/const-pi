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
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fac.h"

typedef struct {
    unsigned fac;
    unsigned pow;
    unsigned nxt;
} sieve_t;

sieve_t  *sieve = NULL;
unsigned sieve_size;
fac_t    ftmp, fmult;

#define INIT_FACS 32

static void fac_reset(fac_t f)
{
    f[0].num_facs = 0;
}

static void fac_resize(fac_t f, unsigned s)
{
    if (f[0].max_facs < s) {
        fac_clear(f);
        fac_init_size(f, s);
    }
}

// remove factors of power 0
static void fac_compact(fac_t f)
{
    unsigned i, j;
    for (i=0, j=0; i<f[0].num_facs; i++) {
        if (f[0].pow[i]>0) {
            if (j<i) {
                f[0].fac[j] = f[0].fac[i];
                f[0].pow[j] = f[0].pow[i];
            }
            j++;
        }
    }
    f[0].num_facs = j;
}

static void fac_mul2(fac_t r, fac_t f, fac_t g)
{
    unsigned i, j, k;

    for (i=j=k=0; i<f[0].num_facs && j<g[0].num_facs; k++) {
        if (f[0].fac[i] == g[0].fac[j]) {
            r[0].fac[k] = f[0].fac[i];
            r[0].pow[k] = f[0].pow[i] + g[0].pow[j];
            i++;
            j++;
        } else if (f[0].fac[i] < g[0].fac[j]) {
            r[0].fac[k] = f[0].fac[i];
            r[0].pow[k] = f[0].pow[i];
            i++;
        } else {
            r[0].fac[k] = g[0].fac[j];
            r[0].pow[k] = g[0].pow[j];
            j++;
        }
    }
    for (; i<f[0].num_facs; i++, k++) {
        r[0].fac[k] = f[0].fac[i];
        r[0].pow[k] = f[0].pow[i];
    }
    for (; j<g[0].num_facs; j++, k++) {
        r[0].fac[k] = g[0].fac[j];
        r[0].pow[k] = g[0].pow[j];
    }
    r[0].num_facs = k;
    assert(k<=r[0].max_facs);
}

////////////////////////////////////////////////////

void fac_init_size(fac_t f, unsigned s)
{
    if (s<INIT_FACS)
        s=INIT_FACS;

    f[0].fac  = (unsigned *)malloc(s*sizeof(unsigned)*2);
    f[0].pow  = f[0].fac + s;
    f[0].max_facs = s;

    fac_reset(f);
}

void fac_init(fac_t f)
{
    fac_init_size(f, INIT_FACS);
}

void fac_clear(fac_t f)
{
    free(f[0].fac);
}

// f = base^pow
void fac_set_bp(fac_t f, unsigned base, unsigned pow)
{
    unsigned i;
    if (!(base < sieve_size)) {
        printf("base %d sieve_size %d\n", base, sieve_size);
        exit(0);
    }
    for (i=0, base/=2; base>0; i++, base = sieve[base].nxt) {
        f[0].fac[i] = sieve[base].fac;
        f[0].pow[i] = sieve[base].pow*pow;
    }
    f[0].num_facs = i;
    assert(i <= f[0].max_facs);
}

// f *= base^pow
void fac_mul_bp(fac_t f, unsigned base, unsigned pow)
{
    fac_set_bp(ftmp, base, pow);
    fac_mul(f, ftmp);
}

// r = f*g
// f *= g
void fac_mul(fac_t f, fac_t g)
{
    fac_t t;
    fac_resize(fmult, f[0].num_facs + g[0].num_facs);
    fac_mul2(fmult, f, g);
    t[0]    = f[0];
    f[0]    = fmult[0];
    fmult[0] = t[0];
}

void fac_show(fac_t f)
{
    unsigned i;
    for (i=0; i<f[0].num_facs; i++)
        if (f[0].pow[i]==1)
            printf("%d ", f[0].fac[i]);
        else
            printf("%d^%d ", f[0].fac[i], f[0].pow[i]);
    printf("\n");
}

/////////////////////////////////////////////////////////////

// convert factorized form to number
static void bs_mul(mpz_t r, fac_t f, unsigned a, unsigned b)
{
    unsigned i, j, m;
    if (b-a<=32) {
        mpz_set_ui(r, 1);
        for (i=a; i<b; i++)
            for (j=0; j<f[0].pow[i]; j++)
                mpz_mul_ui(r, r, f[0].fac[i]);
    } else {
        mpz_t r2;
        mpz_init(r2);
        m = (a + b)/2;
        bs_mul(r2, f, a, m);
        bs_mul(r, f, m, b);
        mpz_mul(r, r, r2);
        mpz_clear(r2);
    }
}

// f /= gcd(f,g), g /= gcd(f,g)
void fac_gcd_compact(fac_t gcd, fac_t f, fac_t g)
{
    unsigned i, j, k, c;

    fac_resize(gcd, min(f[0].num_facs, g[0].num_facs));
    for (i=j=k=0; i<f[0].num_facs && j<g[0].num_facs; ) {
        if (f[0].fac[i] == g[0].fac[j]) {
            c = min(f[0].pow[i], g[0].pow[j]);
            f[0].pow[i] -= c;
            g[0].pow[j] -= c;
            gcd[0].fac[k] = f[0].fac[i];
            gcd[0].pow[k] = c;
            i++;
            j++;
            k++;
        } else if (f[0].fac[i] < g[0].fac[j]) {
            i++;
        } else {
            j++;
        }
    }
    gcd[0].num_facs = k;
    assert(k <= gcd[0].max_facs);
    if (k) {
        fac_compact(f);
        fac_compact(g);
    }
}

void fac_to_mpz(mpz_t r, fac_t f)
{
    bs_mul(r, f, 0, f[0].num_facs);
}


void fac_sieve_init(unsigned size)
{
    unsigned m, n, i, j, k;
    sieve_t  *s;

    if (size == 0)
        return;

    fac_init(ftmp);
    fac_init(fmult);

    sieve = (sieve_t *)malloc(sizeof(sieve_t)*size/2);

    s = sieve;
    n = sieve_size = size;
    m = (unsigned)sqrt(n);
    memset(s, 0, sizeof(sieve_t)*n/2);

    s[1/2].fac = 1;
    s[1/2].pow = 1;

    for (i=3; i<=n; i+=2) {
        if (s[i/2].fac == 0) {
            s[i/2].fac = i;
            s[i/2].pow = 1;
            if (i<=m) {
                for (j=i*i, k=i/2; j<=n; j+=i+i, k++) {
                    if (s[j/2].fac==0) {
                        s[j/2].fac = i;
                        if (s[k].fac == i) {
                            s[j/2].pow = s[k].pow + 1;
                            s[j/2].nxt = s[k].nxt;
                        } else {
                            s[j/2].pow = 1;
                            s[j/2].nxt = k;
                        }
                    }
                }
            }
        }
    }
}

void fac_sieve_clear()
{
    if (sieve)
        free(sieve);
}

