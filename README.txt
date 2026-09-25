PROGRAM FILES SUPPLEMENTING THE PAPER:
An isoperimetric problem for Fourier zeros of centrally symmetric convex bodies

BY: Javier Gómez-Serrano, Michael Levitin, Daniel Platt and Iosif Polterovich


REQUIREMENTS:

    - C++17 compiler and Make
    - FLINT 3 with Arb, GMP and MPFR
    - A platform with 64-bit unsigned long

    The Makefile expects flint/arb.h and links with -lflint -lmpfr -lgmp.


CONTENTS:

    - certify.cc -- Main program to check the numerical bounds and root brackets
    - functions.cc / functions.h -- Formulas for the polygon profile and bounds
    - methods.cc / methods.h -- Interval subdivision and bound checks
    - Makefile -- Builds and runs the certificate

    Formula explanations and references to the paper are in the code comments.


COMPILATION AND RUNNING:

    From the project directory:
        make certify

    To build without running, use 'make'.
    To run with higher precision:
        ./build/certify --bits 512

    The default precision is 384 bits; the minimum is 256 bits.
    Add --quiet to suppress success messages.

    A successful run prints CERTIFICATE PASSED and returns exit status zero.
    Failed or inconclusive checks return a nonzero status.
