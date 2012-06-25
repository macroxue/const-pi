#
# Calculate Mathematical Constants Using GMP (GNU Multiple Precision)
#
# Copyright (C) 2006-2012 Hanhong Xue (macroxue at yahoo dot com)
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA 02111-1307 USA
#

#
# OS dependent options
#
OS= $(shell uname -o)
ifeq ($(OS), Cygwin)
OS_OPT= -Wl,--stack,16777216
endif

#
# Hardware dependent options
#
HARDWARE= $(shell uname -m)
ifeq ($(HARDWARE), i686)
HARDWARE_OPT= -mtune=pentium4 -mno-sse2
else 
ifeq ($(HARDWARE), x86_64)
HARDWARE_OPT= -m64 -mtune=amdfam10
else 
ifeq ($(HARDWARE), ppc64)
HARDWARE_OPT= -mtune=power5
else
HARDWARE_OPT= 
endif
endif
endif

OPT= -O3 $(OS_OPT) $(HARDWARE_OPT) -Wall -static
LIB= -L/usr/local/lib -lgmp -lm
SRC= const.cpp Args.h my.c my.h fac.c fac.h
EXE= const
all: $(EXE)

const: $(SRC)
	g++ $(OPT) $(filter-out %.h,$^) -o $@ $(LIB)

clean:
	rm -f $(EXE)

style:
	astyle -A3 $(SRC)
