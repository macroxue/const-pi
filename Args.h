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

#ifndef _ARGS_H
#define _ARGS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

class Args
{
public:
    Args(int argc, char *argv[]);

    const char *  GetName() const {
        return name;
    }

protected:
    void usage();

protected:
    const char *  name;                 // name of the constant
    size_t        digits;               // number of digits to calculate
    const char *  file;                 // file to save the result
    size_t        gcd_level;            // GCD removal until this level
    size_t        verbose_level;        // verbose level

    static const size_t  MAX_VERBOSE_LEVEL = 10;
};

Args::Args(int argc, char *argv[])
    : name(NULL), digits(1000000), file(NULL), gcd_level(3), verbose_level(0)

{
    int opt;
    while ((opt = getopt(argc, argv, "d:o:g:v:h")) != -1) {
        switch (opt) {
        case 'd':
            digits = atol(optarg);
            switch (optarg[strlen(optarg)-1]) {
            case 'K':
                digits <<= 10;
                break;
            case 'M':
                digits <<= 20;
                break;
            case 'G':
                digits <<= 30;
                break;
            case 'k':
                digits *= 1000;
                break;
            case 'm':
                digits *= 1000000;
                break;
            case 'g':
                digits *= 1000000000;
                break;
            }
            break;
        case 'o':
            file = optarg;
            break;
        case 'g':
            gcd_level = atol(optarg);
            break;
        case 'v':
            verbose_level = atol(optarg);
            break;
        default :
            usage();
            break;
        }
    }

    int pos, i;
    for (pos = 0, i = optind; i < argc; pos++, i++) {
        switch (pos) {
        case 0:
            name = argv[i];
            break;
        default:
            printf("Error: Extra argument '%s'.\n", argv[i]);
            usage();
        }
    }

    if (name == NULL) {
        printf("Specify the name of the constant to compute.\n");
        usage();
    }
    if (digits <= 0) {
        printf("The number of digits must be positive.\n");
        usage();
    }
    if (verbose_level < 0) {
        printf("Verbose level must be positive.\n");
        usage();
    }
    if (verbose_level > MAX_VERBOSE_LEVEL)
        verbose_level = MAX_VERBOSE_LEVEL;
}

void Args::usage()
{
    const char *usage_format =
        "Usage: const name [-d digits] [-o file] [-g level] [-v level]\n"
        "  name        One of the following constants.\n"
        "                pi      : pi using Chudnovsky formula\n"
        "                pi-rama : pi using Ramanujan formula\n"
        "                pi-agm  : pi using AGM\n"
        "                e       : e\n"
        "  -d digits   Number of digits to calculate. Default: %lu.\n"
        "              'K','M','G','k','m','g' notations can be used.\n"
        "  -o file     File to save result. Default: result not saved.\n"
        "  -g level    Remove GCD until this level. Default: %lu.\n"
        "  -v level    Verbose level. Default: %lu.\n";
    printf(usage_format, digits, gcd_level, verbose_level);
    exit(-1);
}

#endif
