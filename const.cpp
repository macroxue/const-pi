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

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include "Args.h"
#include "my.h"
#include "fac.h"

#define p1  (pstack[top])
#define q1  (qstack[top])
#define g1  (gstack[top])
#define fp1 (fpstack[top])
#define fg1 (fgstack[top])
#define p2  (pstack[top+1])
#define q2  (qstack[top+1])
#define g2  (gstack[top+1])
#define fp2 (fpstack[top+1])
#define fg2 (fgstack[top+1])

#define VERBOSE(statement) if (level < verbose_level) statement

const double PI    = 3.14159265358979323846;
const double E     = 2.71828182845904523534;
const double logPI = log10(PI);
const double logE  = log10(E);
const double log_2 = log10(2);

class Constant : public Args
{
protected:
    double log10fac(double n) {
        //log n! ~= (log PI+log 2)/2 + (n+0.5) log n - (n-1/12/n) log E
        return (logPI + log_2)/2 + (n+0.5)*log10(n) - (n-1/12/n)*logE;
    }

protected:
    static const int ALIGN = 20, TAB = 2;
    FILE  *out_fp;
    int   out_count;

    void output(const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        vprintf(fmt, ap);
        fflush(stdout);
        if (out_fp) {
            va_start(ap, fmt);
            vfprintf(out_fp, fmt, ap);
        }
    }

    void output_progress(const char *str) {
        printf("%s", str);
        fflush(stdout);
        out_count += strlen(str);
    }

    void output_retract() {
        for (int i = 0; i < out_count; i++)
            printf("\b \b");
        out_count = 0;
    }

    void output_result(mpf_t result) {
        if (out_fp) {
            output_progress("output to file ...");
            time_begin();
            fputc('\n', out_fp);
            my_out_str(out_fp, 10, digits+2, result);
            fputc('\n', out_fp);
            output_retract();
            time_end("output");
        }
    }

protected:
    clock_t begin[MAX_VERBOSE_LEVEL*2], end[MAX_VERBOSE_LEVEL*2];
    int time_depth;

    void time_begin(const char *step_name = NULL) {
        time_depth++;
        if (step_name)
            output("%*s : ", ALIGN + time_depth * TAB, step_name);
        begin[time_depth] = clock();
    }

    void time_stamp(const char *step_name) {
        end[time_depth] = clock();
        output("%6.3f s\n", time_diff(begin[time_depth], end[time_depth]));
        begin[time_depth] = end[time_depth];
        output("%*s : ", ALIGN + time_depth * TAB, step_name);
    }

    void time_end(const char *step_name = NULL) {
        end[time_depth] = clock();
        if (step_name)
            output("%*s : ", ALIGN + time_depth * TAB, step_name);
        output("%6.3f s\n", time_diff(begin[time_depth], end[time_depth]));
        time_depth--;
    }

    double time_diff(clock_t begin, clock_t end) {
        return (double)(end - begin)/CLOCKS_PER_SEC;
    }

protected:
    long max_rss;

    void check_resource_usage() {
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        if (max_rss < usage.ru_maxrss)
            max_rss = usage.ru_maxrss;
    }

    void show_resource_usage() {
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        output("\n%*s : %.3f s\n", ALIGN, "user time",
               usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1e-6);
        output("%*s : %.3f s\n", ALIGN, "sys time",
               usage.ru_stime.tv_sec + usage.ru_stime.tv_usec*1e-6);
        output("%*s : %ld KB\n", ALIGN, "memory usage", max_rss);
    }

public:
    Constant(Args &args)
        : Args(args), out_fp(NULL), out_count(0), time_depth(-1), max_rss(0) {
        if (file) {
            out_fp = fopen(file, "wt");
            if (!out_fp) {
                printf("Error: cannot open %s for output.", file);
                exit(-1);
            }
        }
        time_t start_time = time(NULL);
        output("Calculation started at %s\n", ctime(&start_time));
    }

    ~Constant() {
        check_resource_usage();
        show_resource_usage();
        time_t finish_time = time(NULL);
        output("\nCalculation finished at %s", ctime(&finish_time));
        if (out_fp)
            fclose(out_fp);
    }
};

class ConstSeries : public Constant
{
protected:
    size_t   terms;
    size_t   depth;
    size_t   top;
    double   progress, percent;
    size_t   progress_count;
    mpz_t    *pstack, *qstack, *gstack;
    fac_t    *fpstack, *fgstack;
    mpz_t    gcd;
    fac_t    fgcd;
    mpf_t    pf, qf;

protected:
    virtual size_t digits_to_terms(size_t d) = 0;
    virtual void  compute_term(size_t b) = 0;
    virtual void  combine_series(size_t a, size_t b, size_t level) = 0;

protected:
    virtual size_t get_sieve_size(size_t terms) {
        return 0;
    }
    virtual double split_ratio(size_t a, size_t b, size_t level) {
        return 0.5;
    }
    virtual void integer_processing() {}
    virtual void floating_point_processing() {}

private:
    size_t get_depth(size_t a, size_t b, size_t level, size_t right) {
        size_t mid;
        if (b-a == 1)
            return level + right;
        else {
            mid = (size_t)(a+(b-a)*split_ratio(a, b, level));
            if (mid - a >= b - mid)
                return get_depth(a, mid, level+1, 0);
            else
                return get_depth(mid, b, level+1, 1);
        }
    }

private:
    void split_series(size_t a, size_t b, size_t level) {
        size_t mid;
        if (b-a == 1) {
            compute_term(b);
            if (b > (size_t)progress) {
                output_progress(".");
                progress += percent*2;
            }
        } else {
            mid = (size_t)(a+(b-a)*split_ratio(a, b, level));
            split_series(a, mid, level+1);
            top++;
            split_series(mid, b, level+1);
            top--;
            if (level < verbose_level) {
                output_retract();
                time_begin();
            }
            combine_series(a, b, level);
            if (level < verbose_level) {
                output("%*s : P[%ld] Q[%ld] G[%ld]\n",
                       ALIGN + (time_depth+1) * TAB, "sizes",
                       mpz_sizeinbase(p1, 10),
                       mpz_sizeinbase(q1, 10), mpz_sizeinbase(g1, 10));
                char level_str[32];
                sprintf(level_str, "level %ld", level);
                time_end(level_str);
            }
        }
    }

public:
    void Calculate() {
        terms = digits_to_terms(digits);
        percent = terms/100.0;

        depth = get_depth(0, terms, 0, 0) + 1;
        output("%ld digits of %s, %ld terms.\n", digits, name, terms);
        time_begin();

        /* initialize stack */
        time_begin("initialization");
        fac_sieve_init( get_sieve_size(terms) );
        pstack = (mpz_t *)malloc(sizeof(mpz_t)*depth);
        qstack = (mpz_t *)malloc(sizeof(mpz_t)*depth);
        gstack = (mpz_t *)malloc(sizeof(mpz_t)*depth);
        fpstack = (fac_t *)malloc(sizeof(fac_t)*depth);
        fgstack = (fac_t *)malloc(sizeof(fac_t)*depth);
        for (size_t i=0; i<depth; i++) {
            mpz_init(pstack[i]);
            mpz_init(qstack[i]);
            mpz_init(gstack[i]);
            fac_init(fpstack[i]);
            fac_init(fgstack[i]);
        }
        mpz_init(gcd);
        fac_init(fgcd);
        my_init(mpf_get_default_prec());
        time_end();

        /* binary splitting */
        time_begin();
        split_series(0, terms, 0);
        output_retract();
        time_end("binary splitting");

        /* free some resources */
        check_resource_usage();
        for (size_t i=1; i<depth; i++) {
            mpz_clear(pstack[i]);
            mpz_clear(qstack[i]);
            mpz_clear(gstack[i]);
            fac_clear(fpstack[i]);
            fac_clear(fgstack[i]);
        }
        mpz_clear(gcd);
        fac_clear(fgcd);
        fac_sieve_clear();

        /* post processing for big integers */
        integer_processing();

        /* prepare to convert integers to floats */
        //mpf_set_default_prec((size_t)(digits*BITS_PER_DIGIT+16));
        mpf_init(pf);
        mpf_init(qf);
        mpf_set_z(pf, p1);
        mpf_set_z(qf, q1);
        mpf_div_2exp(pf, pf, mpz_sizeinbase(q1,2));
        mpf_div_2exp(qf, qf, mpz_sizeinbase(q1,2));

        /* free stacks */
        check_resource_usage();
        mpz_clear(pstack[0]);
        mpz_clear(qstack[0]);
        mpz_clear(gstack[0]);
        fac_clear(fpstack[0]);
        fac_clear(fgstack[0]);
        free(pstack);
        free(qstack);
        free(gstack);
        free(fpstack);
        free(fgstack);

        /* post processing for big floating points */
        //my_init(mpf_get_default_prec());
        floating_point_processing();
        check_resource_usage();
        my_clear();

        /* overall time */
        output("-----------------------------------\n");
        time_end("calculation");

        output_result(pf);
        check_resource_usage();
        mpf_clear(pf);
        mpf_clear(qf);
    }

public:
    ConstSeries(Args &args)
        : Constant(args), top(0), progress(0) {
        mpf_set_default_prec((size_t)(digits*BITS_PER_DIGIT+16));
    }
};

class ConstPi : public ConstSeries
{
protected:
    size_t A, B, C, D, R;

protected:
    void combine_series(size_t a, size_t b, size_t level) {
        /*
           p(a,b) = p(a,m) * p(m,b)
           g(a,b) = g(a,m) * g(m,b)
           q(a,b) = q(a,m) * p(m,b) + q(m,b) * g(a,m)
         */
        if (level >= gcd_level) {
            VERBOSE( time_begin() );
            fac_gcd_compact(fgcd, fp2, fg1);
            if (fgcd[0].num_facs) {
                VERBOSE( time_begin("gcd") );
                fac_to_mpz(gcd, fgcd);
                VERBOSE( time_stamp("p2=p2/gcd") );
                mpz_divexact(p2, p2, gcd);
                VERBOSE( time_stamp("g1=g1/gcd") );
                mpz_divexact(g1, g1, gcd);
                VERBOSE( time_end() );

                VERBOSE(
                    output("%*s : p2[%ld] g1[%ld] gcd[%ld]\n",
                           ALIGN + (time_depth+1) * TAB, "sizes",
                           mpz_sizeinbase(p2, 10),
                           mpz_sizeinbase(g1, 10), mpz_sizeinbase(gcd, 10))
                );
            }
            VERBOSE( time_end("remove gcd") );
        }

        VERBOSE( time_begin("P=p1*p2") );
        mpz_mul(p1, p1, p2);
        VERBOSE( time_stamp("Q=q1*p2") );
        mpz_mul(q1, q1, p2);
        VERBOSE( time_stamp("Q+=q2*g1") );
        mpz_addmul(q1, q2, g1);

        if (b < terms) {
            VERBOSE( time_stamp("G=g1*g2") );
            mpz_mul(g1, g1, g2);
        }
        VERBOSE( time_end() );

        if (level >= gcd_level) {
            fac_mul(fp1, fp2);
            if (b < terms)
                fac_mul(fg1, fg2);
        }
    }

protected:
    void floating_point_processing() {
        time_begin("division");
        mpf_div(qf, pf, qf);
        time_end();

#ifdef INV_SQRT
        time_begin("inverse square root");
        my_invsqrt_ui(pf, R);
#else
        time_begin("square root");
        my_sqrt_ui(pf, R);
#endif
        time_end();

        time_begin("multiplication");
        mpf_mul(pf, qf, pf);
        time_end();
    }

public:
    ConstPi(Args &args)
        : ConstSeries(args) {
    }
};

class ConstPiRama : public ConstPi
{
    /*
            sqrt(8)  inf  (4n)!*(1103 + 26390 n)
    1/pi = --------- Sum ------------------------.
             9801    n=0   (n!)^4 * 396 ^ (4n)
    */
    size_t digits_to_terms(size_t d) {
        /*
          (4n)!(A+Bn)  sqrt(2048)
          ----------- ----------- < 10^(-d)
             (n!)^4    C^(4n+2)

          4log n! + (4n+2)log C - log(4n)! - log(A+Bn) - 0.5*log 2048 > d

          4(n+0.5)log C > d + log(4n)! + log(A+Bn) + 0.5*log 2048 - 4log n!
        */
        double n, n1 = d;
        do {
            n = n1;
            n1 = (d + log10fac(4*n) + log10(A+B*n) + log10(D)/2
                  - 4*log10fac(n)) / log10(C) / 4 - 0.5;
        } while (fabs(n - n1) > 0.01);
        return (size_t)(n1 + 1);
    }

    size_t get_sieve_size(size_t terms) {
        return max(99+1, terms*4);
    }

    double split_ratio(size_t a, size_t b, size_t level) {
        return a == 0 ? 0.5345 : 0.5005;
    }

    void compute_term(size_t b) {
        /*
           g(b-1,b) = (4b-3)(2b-1)(4b-1)
           p(b-1,b) = b^3 * C^4/8
           q(b-1,b) = g(b-1,b)*(A+Bb).
         */
        mpz_set_ui(p1, b);
        mpz_mul_ui(p1, p1, b);
        mpz_mul_ui(p1, p1, b);
        mpz_mul_ui(p1, p1, C/2*C/2*C/2*C);

        mpz_set_ui(g1, 2*b-1);
        mpz_mul_ui(g1, g1, 4*b-1);
        mpz_mul_ui(g1, g1, 4*b-3);

        mpz_set_ui(q1, b);
        mpz_mul_ui(q1, q1, B);
        mpz_add_ui(q1, q1, A);
        mpz_mul   (q1, q1, g1);

        size_t i=b;
        while ((i&1)==0) i>>=1;
        fac_set_bp(fp1, i, 3);    // b^3
        fac_mul_bp(fp1, 99, 4);

        fac_set_bp(fg1, 2*b-1, 1);   // 2b-1
        fac_mul_bp(fg1, 4*b-1, 1);   // 4b-1
        fac_mul_bp(fg1, 4*b-3, 1);   // 4b-3
    }

    void integer_processing() {
        /*
          1/pi = sqrt(2048)*(A+Q/P)/(396^2)

                396^2               (396^2)*P*sqrt(2048)
          pi = ------------------ = --------------------
               sqrt(2048) (A+Q/P)     2048*(A*P+Q)
        */
        mpz_addmul_ui(q1, p1, A);
#ifdef INV_SQRT
#else
        mpz_mul_ui(q1, q1, D);
#endif
        mpz_mul_ui(p1, p1, C*C);
    }

public:
    ConstPiRama(Args &args)
        : ConstPi(args) {
        name = "pi using Ramanujan formula";
        A = 1103;
        B = 26390;
        C = 396;
        D = R = 2048;
    }
};

class ConstPiChud : public ConstPi
{
    /*
                12         inf  (-1)^n * (6n)! * (13591409 + 545140134 n)
    1/pi = --------------- Sum ------------------------------------------.
             640320^(3/2)  n=0  (3n)! * (n!)^3 * 640320^(3n)
    */
    size_t digits_to_terms(size_t d) {
        /*
           (6n)!(A+Bn)   12
           ----------- ---------- < 10^(-d)
           (3n)!(n!)^3 C^(3n+1.5)

           log(3n)! + 3log n! + (3n+1.5)log C - log(6n)! - log(A+Bn) - log12 > d

           3(n+0.5)log C > d + log(6n)! + log(A+Bn) + log12 - log(3n)! - 3log n!
         */
        double n, n1 = d;
        do {
            n = n1;
            n1 = (d + log10fac(6*n) + log10(A+B*n) + log10(D)
                  - log10fac(3*n) - 3*log10fac(n)) / log10(C) / 3 - 0.5;
        } while (fabs(n - n1) > 0.01);
        return (size_t)(n1 + 1);
    }

    size_t get_sieve_size(size_t terms) {
        return max(3*5*23*29+1, terms*6);
    }

    double split_ratio(size_t a, size_t b, size_t level) {
        return a == 0 ? 0.5224 : 0.5048;
    }

    void compute_term(size_t b) {
        /*
           g(b-1,b) = (6b-5)(2b-1)(6b-1)
           p(b-1,b) = b^3 * C^3 / 24
           q(b-1,b) = (-1)^b*g(b-1,b)*(A+Bb).
         */
        mpz_set_ui(p1, b);
        mpz_mul_ui(p1, p1, b);
        mpz_mul_ui(p1, p1, b);
        mpz_mul_ui(p1, p1, (C/24)*(C/24));
        mpz_mul_ui(p1, p1, C*24);

        mpz_set_ui(g1, 2*b-1);
        mpz_mul_ui(g1, g1, 6*b-1);
        mpz_mul_ui(g1, g1, 6*b-5);

        mpz_set_ui(q1, b);
        mpz_mul_ui(q1, q1, B);
        mpz_add_ui(q1, q1, A);
        mpz_mul   (q1, q1, g1);
        if (b%2)
            mpz_neg(q1, q1);

        size_t i=b;
        while ((i&1)==0) i>>=1;
        fac_set_bp(fp1, i, 3);    // b^3
        fac_mul_bp(fp1, 3*5*23*29, 3);
        fp1[0].pow[0]--;

        fac_set_bp(fg1, 2*b-1, 1);   // 2b-1
        fac_mul_bp(fg1, 6*b-1, 1);   // 6b-1
        fac_mul_bp(fg1, 6*b-5, 1);   // 6b-5
    }

    void integer_processing() {
        /*
          1/pi = 12*(A+Q/P)/640320/sqrt(640320)

                   (640320^2)*P               640320*P*sqrt(640320)
          pi = --------------------------- = -------------------------
                 sqrt(640320)*12*(A*P+Q)           12*(A*P+Q)
        */
        mpz_addmul_ui(q1, p1, A);
#ifdef INV_SQRT
        mpz_mul_ui(p1, p1, (C/D)*(C/64));
        mpz_mul_2exp(p1, p1, 6);
#else
        mpz_mul_ui(p1, p1, C/D);
#endif
    }

public:
    ConstPiChud(Args &args)
        : ConstPi(args) {
        name = "pi using Chudnovsky formula";
        A = 13591409;
        B = 545140134;
        C = R = 640320;
        D = 12;
    }
};

class ConstPiAgm : public Constant
{
public:
    ConstPiAgm(Args &args)
        : Constant(args) {
        name = "pi using AGM";
        mpf_set_default_prec((size_t)(digits*BITS_PER_DIGIT+16));
    }

    void Calculate() {
        output("%ld digits of %s.\n", digits, name);
        time_begin();

        /* initialization */
        time_begin("initialization");
        mpf_t a, b, a2, b2, c2, sum;
        mpf_init(a);
        mpf_init(b);
        mpf_init(a2);
        mpf_init(b2);
        mpf_init(c2);
        mpf_init(sum);

        mpf_set_ui(a, 1);
        mpf_set_ui(a2, 1);
        mpf_set_d(b2, 0.5);
        mpf_set_d(c2, 0.5);
        mpf_set_ui(sum, 1);

        long prec0 = mpf_get_default_prec();
        my_init(prec0);

        time_end();

        /* iterations */
        long count, prec;
        for (count = 0, prec = -1; prec < prec0 * 2 + 10;
                count++, prec = prec * 2 + 10) {
            time_begin();

            /*   c2   = a2 - b2; */
            mpf_sub(c2, a2, b2);

            /*   sum -=(1<<n)*c2; */
            mpf_mul_2exp(c2, c2, count);
            mpf_sub(sum, sum, c2);

            /*  b    = sqrt(b2); */
            if (verbose_level) time_begin("square root");
            my_sqrt(b, b2);
            if (verbose_level) time_end();

            /*  a    =(a+b)/2; */
            mpf_add(a, a, b);
            mpf_div_2exp(a, a, 1);

            /*  c2   =(a-b)*(a-b); */
            mpf_sub(c2, a, b);
            if (verbose_level) time_begin("square");
            mpf_mul(c2, c2, c2);
            if (verbose_level) time_end();

            /*  b2   =((a2+b2)/4-c2)*2; */
            mpf_add(b2, b2, a2);
            mpf_div_2exp(b2, b2, 2);
            mpf_sub(b2, b2, c2);
            mpf_mul_2exp(b2, b2, 1);

            /*  a2   = b2+c2; */
            mpf_add(a2, b2, c2);

            char message[50];
            sprintf(message, "iteration %ld", count+1);
            time_end(message);
        }
        time_begin("division");
        mpf_mul_2exp(a2, a2, 1);
        mpf_div(a2, a2, sum);
        time_end();

        output("-----------------------------------\n");
        time_end("calculation");
        output_result(a2);
    }
};

class ConstE : public ConstSeries
{
    size_t digits_to_terms(size_t d) {
        /*
                 ________               1       1
          n! ~= V 2 PI n  n^n exp(-n + --- + O(---))
                                       12n     n^2

          log n! ~= (log PI+log 2)/2 + (n+0.5) log n - (n-1/12/n) log E

          1/n! < 10^(-d) ==> n! > 10^d ==> log n! > d

        */
        double n, n1 = d;
        do {
            n = n1;
            n1 = (d - (logPI + log_2)/2 + (n-1/12/n)*logE) / log10(n) - 0.5;
        } while (fabs(n - n1) > 0.01);
        return (size_t)(n1 + 1);
    }

    void compute_term(size_t b) {
        mpz_set_ui(p1, b);
        mpz_set_ui(q1, 1);
    }

    void combine_series(size_t a, size_t b, size_t level) {
        VERBOSE( time_begin("P=p1*p2") );
        mpz_mul(p1, p1, p2);
        VERBOSE( time_stamp("Q=q1*p2+q2") );
        mpz_mul(q1, q1, p2);
        mpz_add(q1, q1, q2);
        VERBOSE( time_end() );
    }

    void floating_point_processing() {
        time_begin("division");
        mpf_div(pf, qf, pf);
        mpf_add_ui(pf, pf, 1);
        time_end();
    }

public:
    ConstE(Args &args)
        : ConstSeries(args) {
        name = "e";
    }
};

int main(int argc, char *argv[])
{
    Args args(argc, argv);

    const char *name = args.GetName();
    if (strcasecmp(name, "pi")==0) {
        ConstPiChud pi(args);
        pi.Calculate();
    } else if (strcasecmp(name, "pi-rama")==0) {
        ConstPiRama pi(args);
        pi.Calculate();
    } else if (strcasecmp(name, "pi-agm")==0) {
        ConstPiAgm pi(args);
        pi.Calculate();
    } else if (strcasecmp(name, "e")==0) {
        ConstE e(args);
        e.Calculate();
    } else {
        printf("Unknown constant '%s'\n", name);
        return 1;
    }

    return 0;
}
