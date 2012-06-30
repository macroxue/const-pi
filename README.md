Calculating Mathematical Constants
==================================

This is a program calculating matchematical constants, such as pi and e,
to arbitrary precision using the GNU Multiple Precision (GMP) Arithmetic
Library from http://gmplib.org.

Currently there are two types of algorithms implemented.
 * Series calculation with binary splitting
 * Iterative calculation with Arithmetic-Geometric Mean (AGM) 

Series calculation with binary splitting is generic enough to handle different
series such as Chudnovsky formula for pi, Ramanujan formula for pi and Newton's
series for e.

For now, AGM is only used in calculating pi.


Compiling
---------

Type `make` to compile the program under Linux and other Unix-like operating
systems. The resulting executable is `const`.

If GMP library is not installed in a standard location, such as `/usr/lib`
or `/usr/local/lib`, the `LIB` line in `makefile` needs to be updated to pick
up the new location.


Running
-------

Type `const` without arguments for help and follow the instructions.

