// Implementation file
// Functions needed to certify the polygon zeros and the numerical bounds.

#include "functions.h"

#include "flint/arb_hypgeom.h"

// res = 2*J_1(x)/x = 0F1(;2;-x^2/4).
// The value at zero is 1, so the positivity check can start at x=0.
void arb_normalized_bessel_j1(arb_t res, const arb_t x, const void *params, slong prec){
  (void) params;
  arb_t a, z;
  arb_init(a); arb_init(z);
  arb_set_ui(a,2);
  arb_mul(z,x,x,prec);
  arb_div_ui(z,z,4,prec);
  arb_neg(z,z);
  arb_hypgeom_0f1(res,a,z,0,prec);
  arb_clear(a); arb_clear(z);
}

// Equations (4.2)-(4.3): q_n for the area-pi polygon P_{2n}, n >= 2.
void arb_polygon_q(arb_t res, const arb_t rho, ulong n, slong prec){
  arb_t pi, delta, R, phi, prefactor, apothem, transverse;
  arb_t alpha, cosine, sine, arg1, arg2, term, aux, sum;
  arb_init(pi); arb_init(delta); arb_init(R); arb_init(phi);
  arb_init(prefactor); arb_init(apothem); arb_init(transverse);
  arb_init(alpha); arb_init(cosine); arb_init(sine);
  arb_init(arg1); arb_init(arg2); arb_init(term); arb_init(aux); arb_init(sum);

  arb_const_pi(pi,prec);
  arb_div_ui(delta,pi,2*n,prec);
  // R_n^2 = 2*delta/sin(2*delta).
  arb_mul_ui(aux,delta,2,prec);
  arb_sin(R,aux,prec);
  arb_div(R,aux,R,prec);
  arb_sqrt(R,R,prec);

  // Minimising direction: phi*=0 for even n, phi*=delta for odd n.
  if(n%2 == 0) arb_zero(phi);
  else arb_set(phi,delta);

  // Simplify 4*R_n^2*cos(delta)*sin(delta) exactly to 2*pi/n.
  arb_mul_ui(prefactor,pi,2,prec);
  arb_div_ui(prefactor,prefactor,n,prec);
  arb_cos(aux,delta,prec);
  arb_mul(apothem,R,aux,prec);
  arb_sin(aux,delta,prec);
  arb_mul(transverse,R,aux,prec);
  arb_zero(sum);

  for(ulong k = 0; k < n; k++){
    arb_mul_ui(alpha,pi,k,prec);
    arb_div_ui(alpha,alpha,n,prec);
    arb_sub(alpha,phi,alpha,prec);
    arb_cos(cosine,alpha,prec);
    arb_sin(sine,alpha,prec);

    arb_mul(arg1,rho,apothem,prec);
    arb_mul(arg1,arg1,cosine,prec);
    arb_mul(arg2,rho,transverse,prec);
    arb_mul(arg2,arg2,sine,prec);
    // sinc(x)=sin(x)/x, with sinc(0)=1.
    arb_sinc(arg1,arg1,prec);
    arb_sinc(arg2,arg2,prec);
    arb_mul(term,cosine,cosine,prec);
    arb_mul(term,term,arg1,prec);
    arb_mul(term,term,arg2,prec);
    arb_add(sum,sum,term,prec);
  }
  arb_mul(res,prefactor,sum,prec);

  arb_clear(pi); arb_clear(delta); arb_clear(R); arb_clear(phi);
  arb_clear(prefactor); arb_clear(apothem); arb_clear(transverse);
  arb_clear(alpha); arb_clear(cosine); arb_clear(sine);
  arb_clear(arg1); arb_clear(arg2); arb_clear(term); arb_clear(aux); arb_clear(sum);
}

// Interval callback for q_n; params contains n.
void arb_polygon_q_interval(arb_t res, const arb_t rho, const void *params, slong prec){
  ulong n = *((const ulong *) params);
  arb_polygon_q(res,rho,n,prec);
}

// -D'(rho)=2*pi*J_2(rho)/rho for the disk profile; used away from rho=0.
void arb_J2_quotient(arb_t res, const arb_t rho, const void *params, slong prec){
  (void) params;
  arb_t pi, J2, order;
  arb_init(pi); arb_init(J2); arb_init(order);
  arb_const_pi(pi,prec);
  arb_set_ui(order,2);
  arb_hypgeom_bessel_j(J2,order,rho,prec);
  arb_mul(res,pi,J2,prec);
  arb_mul_ui(res,res,2,prec);
  arb_div(res,res,rho,prec);
  arb_clear(pi); arb_clear(J2); arb_clear(order);
}

// g(s)=s*J_1(j11*s)/j11.
// res = g''''(s) = -2*j11^2*J_0(j11*s)+(j11/s+j11^3*s)*J_1(j11*s).
// params contains j11; the interval for s must exclude zero.
void arb_g_fourth(arb_t res, const arb_t s, const void *params, slong prec){
  arb_srcptr j11 = (arb_srcptr) params;
  arb_t arg, J0, J1, j11_sq, aux, order;
  arb_init(arg); arb_init(J0); arb_init(J1); arb_init(j11_sq); arb_init(aux);
  arb_init(order);

  arb_mul(arg,j11,s,prec);
  arb_zero(order);
  arb_hypgeom_bessel_j(J0,order,arg,prec);
  arb_one(order);
  arb_hypgeom_bessel_j(J1,order,arg,prec);
  arb_mul(j11_sq,j11,j11,prec);

  arb_mul(res,j11_sq,J0,prec);
  arb_mul_si(res,res,-2,prec);
  arb_div(aux,j11,s,prec);
  arb_mul(aux,aux,J1,prec);
  arb_add(res,res,aux,prec);
  arb_mul(aux,j11_sq,j11,prec);
  arb_mul(aux,aux,s,prec);
  arb_mul(aux,aux,J1,prec);
  arb_add(res,res,aux,prec);

  arb_clear(arg); arb_clear(J0); arb_clear(J1); arb_clear(j11_sq); arb_clear(aux);
  arb_clear(order);
}

// gamma(s)=-s^2*J_2(j11*s)/j11.
// res = gamma''(s) = -s*J_1(j11*s)-j11*s^2*J_0(j11*s).
// params contains j11.
void arb_gamma_second(arb_t res, const arb_t s, const void *params, slong prec){
  arb_srcptr j11 = (arb_srcptr) params;
  arb_t arg, J0, J1, aux, order;
  arb_init(arg); arb_init(J0); arb_init(J1); arb_init(aux); arb_init(order);

  arb_mul(arg,j11,s,prec);
  arb_zero(order);
  arb_hypgeom_bessel_j(J0,order,arg,prec);
  arb_one(order);
  arb_hypgeom_bessel_j(J1,order,arg,prec);
  arb_mul(res,s,J1,prec);
  arb_neg(res,res);
  arb_mul(aux,s,s,prec);
  arb_mul(aux,aux,j11,prec);
  arb_mul(aux,aux,J0,prec);
  arb_sub(res,res,aux,prec);

  arb_clear(arg); arb_clear(J0); arb_clear(J1); arb_clear(aux); arb_clear(order);
}

// Jbar_nu(y) = 1-y/(nu+1)+y^2/(2*(nu+1)*(nu+2))
//             -y^3/(6*(nu+1)*(nu+2)*(nu+3)).
// Lower bound for the normalized Bessel function when nu>=6, 0<=y<=pi^2
// (Lemma B.4).
void arb_jbar_polynomial(arb_t res, ulong nu, const arb_t y, slong prec){
  arb_t y2, y3, aux;
  arb_init(y2); arb_init(y3); arb_init(aux);
  arb_mul(y2,y,y,prec);
  arb_mul(y3,y2,y,prec);
  arb_one(res);
  arb_div_ui(aux,y,nu+1,prec); arb_sub(res,res,aux,prec);
  arb_div_ui(aux,y2,2*(nu+1)*(nu+2),prec); arb_add(res,res,aux,prec);
  arb_div_ui(aux,y3,6*(nu+1)*(nu+2)*(nu+3),prec); arb_sub(res,res,aux,prec);
  arb_clear(y2); arb_clear(y3); arb_clear(aux);
}

// Appendix A.1 factors; the six output slots are documented in functions.h.
void arb_master_components(arb_ptr values, ulong n, slong prec){
  arb_t pi, delta, cos_delta, den, R1, R2, aux, y;
  arb_init(pi); arb_init(delta); arb_init(cos_delta); arb_init(den);
  arb_init(R1); arb_init(R2); arb_init(aux); arb_init(y);

  arb_const_pi(pi,prec);
  arb_div_ui(delta,pi,2*n,prec);
  arb_cos(cos_delta,delta,prec);
  // R1=cos(delta)/cos(3*delta/4), R2=cos(delta)/cos(delta/4).
  arb_mul_ui(aux,delta,3,prec); arb_div_ui(aux,aux,4,prec);
  arb_cos(den,aux,prec); arb_div(R1,cos_delta,den,prec);
  arb_div_ui(aux,delta,4,prec);
  arb_cos(den,aux,prec); arb_div(R2,cos_delta,den,prec);

  arb_pow_ui(values,pi,2*n,prec);
  arb_fac_ui(aux,2*n,prec); arb_mul(values,values,aux,prec);
  arb_fac_ui(aux,4*n,prec); arb_div(values,values,aux,prec);

  arb_sqrt_ui(values+1,2,prec);
  arb_mul_ui(values+1,values+1,4*n+2,prec);
  arb_set_str(den,"32.08",prec); arb_div(values+1,values+1,den,prec);
  arb_sub(values+2,R1,R2,prec);
  arb_pow_ui(values+3,R2,2*n+1,prec);
  arb_mul(y,pi,pi,prec); arb_mul(aux,R1,R1,prec); arb_mul(y,y,aux,prec);
  arb_jbar_polynomial(values+4,2*n,y,prec);
  arb_mul(values+5,values+1,values+2,prec);
  arb_mul(values+5,values+5,values+3,prec);
  arb_mul(values+5,values+5,values+4,prec);

  arb_clear(pi); arb_clear(delta); arb_clear(cos_delta); arb_clear(den);
  arb_clear(R1); arb_clear(R2); arb_clear(aux); arb_clear(y);
}
