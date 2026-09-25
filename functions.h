// Header file
// Functions needed to certify the polygon zeros and the numerical bounds.

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "flint/arb.h"

void arb_normalized_bessel_j1(arb_t res, const arb_t x, const void *params, slong prec);
void arb_polygon_q(arb_t res, const arb_t rho, ulong n, slong prec);
void arb_polygon_q_interval(arb_t res, const arb_t rho, const void *params, slong prec);
void arb_J2_quotient(arb_t res, const arb_t rho, const void *params, slong prec);
void arb_g_fourth(arb_t res, const arb_t s, const void *params, slong prec);
void arb_gamma_second(arb_t res, const arb_t s, const void *params, slong prec);
void arb_jbar_polynomial(arb_t res, ulong nu, const arb_t y, slong prec);

// Caller allocates six Arb entries, in the Appendix A.1 table order:
// [0] pi^(2n)*(2n)!/(4n)!; [1] sqrt(2)*(4n+2)/32.08;
// [2] R1-R2; [3] R2^(2n+1); [4] Jbar_{2n}(pi^2*R1^2);
// [5] S_n, the product of entries [1] through [4].
void arb_master_components(arb_ptr values, ulong n, slong prec);

#endif
