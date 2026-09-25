// Instructions: ./build/certify [--bits N] [--quiet]
// Output: PASSED/FAILED for each group of checks; nonzero exit status on failure.
// Checks the numerical bounds in Sections 3, 4 and Appendix A.

#include <algorithm>
#include <iostream>
#include <string.h>
#include <stdlib.h>

#include "flint/arb.h"
#include "flint/arb_hypgeom.h"
#include "flint/flint.h"
#include "flint/fmpq.h"
#include "flint/ulong_extras.h"
#include "functions.h"
#include "methods.h"

using namespace std;

struct root_bracket {
  ulong n;
  const char *lower_exact;
  const char *upper_exact;
};

// Table 1 (tab:enclosures): candidate endpoints of width 10^-12.
// certify_polygon_roots proves these enclosures before later checks use them.
const root_bracket polygon_root_brackets[] = {
  {2,"3544907701811/1000000000000","3544907701812/1000000000000"},
  {3,"3809251227455/1000000000000","3809251227456/1000000000000"},
  {4,"3830839192487/1000000000000","3830839192488/1000000000000"},
  {5,"3831700570935/1000000000000","3831700570936/1000000000000"},
  {6,"3831712243963/1000000000000","3831712243964/1000000000000"},
  {7,"3831708614742/1000000000000","3831708614743/1000000000000"},
  {8,"3831707147744/1000000000000","3831707147745/1000000000000"},
  {9,"3831706546813/1000000000000","3831706546814/1000000000000"},
  {10,"3831706275027/1000000000000","3831706275028/1000000000000"},
  {11,"3831706141601/1000000000000","3831706141602/1000000000000"},
  {12,"3831706071593/1000000000000","3831706071594/1000000000000"},
  {13,"3831706032783/1000000000000","3831706032784/1000000000000"},
  {14,"3831706010248/1000000000000","3831706010249/1000000000000"},
  {15,"3831705996637/1000000000000","3831705996638/1000000000000"},
  {16,"3831705988130/1000000000000","3831705988131/1000000000000"},
  {17,"3831705982652/1000000000000","3831705982653/1000000000000"},
  {18,"3831705979032/1000000000000","3831705979033/1000000000000"},
  {19,"3831705976582/1000000000000","3831705976583/1000000000000"},
  {20,"3831705974891/1000000000000","3831705974892/1000000000000"},
  {21,"3831705973700/1000000000000","3831705973701/1000000000000"}
};

const ulong number_polygon_root_brackets =
  sizeof(polygon_root_brackets)/sizeof(polygon_root_brackets[0]);

const char j11_lower_exact[] = "3831705970207512315/1000000000000000000";
const char j11_upper_exact[] = "3831705970207512316/1000000000000000000";

int print_output = 1;

void passed(const char *name){
  if(print_output) cout << name << ": PASSED" << endl;
}

// Lemma A.1(i)-(ii): bounds on j11, J_0(j11), B0 and C0.
// j11 and C0 are output parameters.
void certify_constants(arb_t j11, arb_t C0, slong prec){
  arb_t left, right, res, bound, pi, aux, J0, B0, order;
  fmpq_t q;
  arb_init(left); arb_init(right); arb_init(res); arb_init(bound);
  arb_init(pi); arb_init(aux); arb_init(J0); arb_init(B0);
  arb_init(order);
  fmpq_init(q);

  fmpq_set_str(q,j11_lower_exact,10);
  fmpq_canonicalise(q);
  arb_set_fmpq(left,q,prec);
  fmpq_set_str(q,j11_upper_exact,10);
  fmpq_canonicalise(q);
  arb_set_fmpq(right,q,prec);
  arb_one(order);
  arb_hypgeom_bessel_j(res,order,left,prec);
  check_positive(res,"J1 at left endpoint");
  arb_hypgeom_bessel_j(res,order,right,prec);
  check_negative(res,"J1 at right endpoint");

  // The sign change gives a zero; positivity before the bracket makes it first.
  subdivision_stats stats;
  if(!test_positivity(arb_normalized_bessel_j1,NULL,"0",j11_lower_exact,
                      120,prec,stats))
    failed("J1 before its first zero");

  arb_union(j11,left,right,prec);
  arb_zero(order);
  arb_hypgeom_bessel_j(J0,order,j11,prec);
  arb_set_str(left,"-0.40275948",prec);
  arb_set_str(right,"-0.40275931",prec);
  check_inside(J0,left,right,"J0(j11)");

  arb_mul(aux,j11,j11,prec);
  arb_abs(res,J0);
  arb_mul(res,res,aux,prec);
  arb_div_ui(res,res,6,prec);
  arb_set_str(bound,"0.9856",prec);
  check_less(res,bound,"j11^2 |J0(j11)|/6");

  // B0 = -j11^2*J_0(j11)*pi^7/90720.
  arb_const_pi(pi,prec);
  arb_pow_ui(res,pi,7,prec);
  arb_mul(B0,aux,J0,prec);
  arb_mul(B0,B0,res,prec);
  arb_neg(B0,B0);
  arb_div_ui(B0,B0,90720,prec);
  arb_set_str(left,"0.196868",prec);
  arb_set_str(right,"0.196869",prec);
  check_inside(B0,left,right,"B0");

  // C0 = j11^3*pi^6/181440, reused in the bridge and effective estimates.
  arb_pow_ui(aux,j11,3,prec);
  arb_pow_ui(res,pi,6,prec);
  arb_mul(C0,aux,res,prec);
  arb_div_ui(C0,C0,181440,prec);
  arb_set_ui(left,149); arb_div_ui(left,left,500,prec);
  arb_set_str(right,"0.2981",prec);
  check_inside(C0,left,right,"C0");

  if(print_output)
    cout << "CERTIFIED CONSTANTS: PASSED (" << stats.boxes
         << " boxes, depth " << stats.depth << ")" << endl;

  arb_clear(left); arb_clear(right); arb_clear(res); arb_clear(bound);
  arb_clear(pi); arb_clear(aux); arb_clear(J0); arb_clear(B0);
  arb_clear(order);
  fmpq_clear(q);
}

// Lemma A.1(iii)-(vi): bounds on J_2, g'''' and gamma''.
void certify_bessel_bounds(const arb_t j11, slong prec){
  arb_t res, bound, aux, pi, order;
  arb_init(res); arb_init(bound); arb_init(aux); arb_init(pi); arb_init(order);

  // Check [j11-0.1,j11+0.1] is contained in [3.731,3.932].
  arb_set_str(bound,"0.1",prec);
  arb_sub(res,j11,bound,prec);
  arb_set_str(bound,"3.731",prec);
  check_greater(res,bound,"left endpoint of Bessel interval");
  arb_set_str(bound,"0.1",prec);
  arb_add(res,j11,bound,prec);
  arb_set_str(bound,"3.932",prec);
  check_less(res,bound,"right endpoint of Bessel interval");

  subdivision_stats stats1, stats2, stats3;
  arb_set_ui(bound,3); arb_div_ui(bound,bound,5,prec);
  if(!test_lower_bound(arb_J2_quotient,NULL,"3731/1000","983/250",bound,
                       120,prec,stats1))
    failed("2*pi*J2(rho)/rho > 3/5");

  arb_set_ui(order,2);
  arb_hypgeom_bessel_j(res,order,j11,prec);
  arb_const_pi(pi,prec);
  arb_mul(res,res,pi,prec);
  arb_mul_ui(res,res,2,prec);
  arb_mul(aux,j11,j11,prec);
  arb_div(res,res,aux,prec);
  arb_set_ui(bound,17); arb_div_ui(bound,bound,100,prec);
  check_greater(res,bound,"2*pi*J2(j11)/j11^2");

  arb_set_ui(bound,13);
  if(!test_absolute_bound(arb_g_fourth,j11,"99/100","101/100",bound,
                          120,prec,stats2))
    failed("fourth derivative of g");
  arb_set_ui(bound,2);
  if(!test_absolute_bound(arb_gamma_second,j11,"99/100","101/100",bound,
                          120,prec,stats3))
    failed("second derivative of gamma");

  if(print_output)
    cout << "BESSEL BOUNDS: PASSED (" << stats1.boxes << ", "
         << stats2.boxes << ", " << stats3.boxes << " boxes)" << endl;

  arb_clear(res); arb_clear(bound); arb_clear(aux); arb_clear(pi); arb_clear(order);
}

// Lemma A.2 and Lemma 4.7: series tails and remainder constants at t=2/25.
// Nonnegative series coefficients extend the tail bounds to smaller arguments.
void certify_elementary_bounds(slong prec){
  arb_t t, secant, secant_tail, cotangent_tail, a, b, c, bound;
  arb_init(t); arb_init(secant); arb_init(secant_tail); arb_init(cotangent_tail);
  arb_init(a); arb_init(b); arb_init(c); arb_init(bound);
  arb_set_str(t,"0.08",prec);

  // secant_tail = sec(t)-1-t^2/2-5*t^4/24. Check 0<secant_tail<t^6/11.
  arb_cos(a,t,prec);
  arb_inv(secant,a,prec);
  arb_mul(a,t,t,prec); arb_div_ui(a,a,2,prec);
  arb_pow_ui(b,t,4,prec); arb_mul_ui(b,b,5,prec); arb_div_ui(b,b,24,prec);
  arb_sub_ui(secant_tail,secant,1,prec);
  arb_sub(secant_tail,secant_tail,a,prec);
  arb_sub(secant_tail,secant_tail,b,prec);
  check_positive(secant_tail,"secant tail");
  arb_pow_ui(a,t,6,prec);
  arb_div(a,secant_tail,a,prec);
  arb_one(bound); arb_div_ui(bound,bound,11,prec);
  check_less(a,bound,"secant tail coefficient");

  // cotangent_tail = 1-t^2/3-t^4/45-t*cot(t). Check 0<tail<t^6/460.
  arb_mul(a,t,t,prec); arb_div_ui(a,a,3,prec);
  arb_pow_ui(b,t,4,prec); arb_div_ui(b,b,45,prec);
  arb_cot(c,t,prec); arb_mul(c,c,t,prec);
  arb_one(cotangent_tail);
  arb_sub(cotangent_tail,cotangent_tail,a,prec);
  arb_sub(cotangent_tail,cotangent_tail,b,prec);
  arb_sub(cotangent_tail,cotangent_tail,c,prec);
  check_positive(cotangent_tail,"cotangent tail");
  arb_pow_ui(a,t,6,prec);
  arb_div(a,cotangent_tail,a,prec);
  arb_one(bound); arb_div_ui(bound,bound,460,prec);
  check_less(a,bound,"cotangent tail coefficient");

  arb_one(a); arb_div_ui(a,a,920,prec);
  arb_one(b); arb_div_ui(b,b,540,prec); arb_add(a,a,b,prec);
  arb_one(b); arb_div_ui(b,b,5000,prec); arb_add(a,a,b,prec);
  arb_one(b); arb_div_ui(b,b,213,prec); arb_add(a,a,b,prec);
  arb_one(bound); arb_div_ui(bound,bound,125,prec);
  check_less(a,bound,"remainder in the square root expansion");

  arb_one(a); arb_div_ui(a,a,124,prec);
  arb_one(b); arb_div_ui(b,b,11,prec); arb_add(a,a,b,prec);
  arb_one(b); arb_div_ui(b,b,24,prec); arb_add(a,a,b,prec);
  arb_one(b); arb_div_ui(b,b,79,prec); arb_add(a,a,b,prec);
  arb_one(bound); arb_div_ui(bound,bound,6,prec);
  check_less(a,bound,"remainder in the profile expansion");

  arb_mul(a,t,t,prec);
  arb_set_str(b,"2.99",prec); arb_div(a,a,b,prec);
  arb_one(bound); arb_div_ui(bound,bound,400,prec);
  check_less(a,bound,"bound for y");

  arb_set_str(a,"2.99",prec); arb_pow_ui(a,a,3,prec);
  arb_mul_ui(a,a,8,prec); arb_inv(a,a,prec);
  arb_one(bound); arb_div_ui(bound,bound,213,prec);
  check_less(a,bound,"cubic bound for y");

  arb_mul(a,t,t,prec); arb_div_ui(a,a,4000,prec);
  arb_one(bound); arb_div_ui(bound,bound,5000,prec);
  check_less(a,bound,"eighth order remainder");

  arb_div_ui(a,secant,125,prec);
  arb_one(bound); arb_div_ui(bound,bound,124,prec);
  check_less(a,bound,"sec(2/25)/125");

  arb_sub_ui(a,secant,1,prec);
  arb_mul(b,t,t,prec); arb_mul_ui(b,b,40,prec);
  arb_div(a,a,b,prec);
  arb_one(bound); arb_div_ui(bound,bound,79,prec);
  check_less(a,bound,"profile secant remainder");

  arb_set_ui(a,1); arb_div_ui(a,a,10,prec);
  arb_set_ui(b,4); arb_div_ui(b,b,625,prec); arb_div_ui(b,b,6,prec);
  arb_add(a,a,b,prec);
  arb_one(bound); arb_div_ui(bound,bound,9,prec);
  check_less(a,bound,"last profile bound");

  passed("ELEMENTARY BOUNDS");
  arb_clear(t); arb_clear(secant); arb_clear(secant_tail); arb_clear(cotangent_tail);
  arb_clear(a); arb_clear(b); arb_clear(c); arb_clear(bound);
}

// Proposition 3.6 and Appendix A.1: check the master inequality.
void certify_master(slong prec){
  // Same column order as arb_master_components: the first entry is rounded up,
  // and all five right-hand entries are rounded down.
  const char *table[3][6] = {
    {"0.0014451", "0.61717", "0.063880", "0.38800", "0.21692", "0.0033181"},
    {"0.000018286", "0.79351", "0.037101", "0.51215", "0.31693", "0.0047785"},
    {"0.00000013969", "0.96984", "0.024083", "0.59568", "0.39758", "0.0055315"}
  };
  const char *ratios[3] = {"2.29", "261", "39000"};
  arb_ptr values = _arb_vec_init(6);
  arb_t bound, pi, lambda, one_minus, a, b, c, d;
  arb_init(bound); arb_init(pi); arb_init(lambda); arb_init(one_minus);
  arb_init(a); arb_init(b); arb_init(c); arb_init(d);

  // The cases n=3,4,5.
  for(ulong n = 3; n <= 5; n++){
    arb_master_components(values,n,prec);
    check_positive(values+4,"positivity of Jbar",n);
    check_less(values,values+5,"master inequality",n);
    arb_set_str(bound,table[n-3][0],prec);
    check_less(values,bound,"left entry in table",n);
    for(int j = 1; j < 6; j++){
      arb_set_str(bound,table[n-3][j],prec);
      check_greater(values+j,bound,"table entry",n);
    }

    // Check the product of the four printed factors and the printed ratio S_n/LHS.
    arb_one(a);
    for(int j = 1; j < 5; j++){
      arb_set_str(b,table[n-3][j],prec);
      arb_mul(a,a,b,prec);
    }
    arb_set_str(bound,table[n-3][5],prec);
    check_greater(a,bound,"product of printed factors",n);
    arb_set_str(b,table[n-3][0],prec);
    arb_div(a,bound,b,prec);
    arb_set_str(bound,ratios[n-3],prec);
    check_greater(a,bound,"ratio of printed bounds",n);
  }

  // Lemma 3.3: lambda=0.0002 bounds the Fourier ratio for n>=3.
  // Check 4/(1-lambda)+4*lambda/(1-lambda)^2+lambda*(1+lambda)/(1-lambda)^3 < 4.01.
  arb_const_pi(pi,prec);
  arb_set_str(lambda,"0.0002",prec);
  arb_div_ui(a,pi,13,prec); arb_pow_ui(a,a,6,prec);
  check_less(a,lambda,"Fourier tail ratio");
  arb_one(one_minus); arb_sub(one_minus,one_minus,lambda,prec);

  arb_set_ui(a,4); arb_div(a,a,one_minus,prec);
  arb_set_str(bound,"4.0009",prec);
  check_less(a,bound,"geometric sum");

  arb_mul_ui(b,lambda,4,prec);
  arb_mul(c,one_minus,one_minus,prec); arb_div(b,b,c,prec);
  arb_set_str(bound,"0.0009",prec);
  check_less(b,bound,"linear geometric sum");

  arb_add_ui(c,lambda,1,prec); arb_mul(c,c,lambda,prec);
  arb_pow_ui(d,one_minus,3,prec); arb_div(c,c,d,prec);
  arb_set_str(bound,"0.0003",prec);
  check_less(c,bound,"quadratic geometric sum");
  arb_add(d,a,b,prec); arb_add(d,d,c,prec);
  arb_set_str(bound,"4.01",prec);
  check_less(d,bound,"total geometric sum");

  // Appendix A.1: check the constants giving S_n > 0.00568/n for n>=6.
  // The paper reduces these bounds to their values at n=6.
  arb_div_ui(a,pi,12,prec); arb_cos(a,a,prec);
  arb_set_str(bound,"0.9659",prec);
  check_greater(a,bound,"cos(pi/12)");

  arb_mul(a,pi,pi,prec); arb_mul_ui(a,a,13,prec); arb_div_ui(a,a,288,prec);
  arb_one(b); arb_sub(a,b,a,prec);
  arb_set_str(bound,"0.5544",prec);
  check_greater(a,bound,"radial power bound");

  arb_mul(a,pi,pi,prec); arb_div_ui(a,a,13,prec);
  arb_one(b); arb_sub(a,b,a,prec);
  arb_set_str(bound,"0.2407",prec);
  check_greater(a,bound,"Jbar lower bound");

  arb_sqrt_ui(a,2,prec); arb_set_str(b,"32.08",prec); arb_div(a,a,b,prec);
  arb_set_str(bound,"0.04408",prec);
  check_greater(a,bound,"constant in master bound");

  arb_set_str(a,"0.04408",prec); arb_set_str(b,"0.9659",prec); arb_mul(a,a,b,prec);
  arb_set_str(b,"0.5544",prec); arb_mul(a,a,b,prec);
  arb_set_str(b,"0.2407",prec); arb_mul(a,a,b,prec);
  arb_set_str(bound,"0.00568",prec);
  check_greater(a,bound,"product in master bound");

  arb_mul(a,pi,pi,prec);
  arb_set_str(bound,"9.8697",prec);
  check_less(a,bound,"pi^2");

  // Lambda_6 < (0.058401)^6 < 3.97e-8.
  arb_set_str(b,"0.058401",prec); arb_pow_ui(b,b,6,prec);
  arb_div_ui(a,a,169,prec); arb_pow_ui(a,a,6,prec);
  check_less(a,b,"Lambda_6 comparison");
  arb_set_str(bound,"0.0000000397",prec);
  check_less(b,bound,"Lambda_6");
  arb_mul_ui(a,bound,6,prec);
  arb_set_str(bound,"0.00000024",prec);
  check_less(a,bound,"W_6");

  arb_mul(a,pi,pi,prec); arb_mul_ui(a,a,7,prec); arb_div_ui(a,a,1350,prec);
  arb_set_str(bound,"0.052",prec);
  check_less(a,bound,"quotient W_(n+1)/W_n");

  passed("MASTER INEQUALITY");
  _arb_vec_clear(values,6);
  arb_clear(bound); arb_clear(pi); arb_clear(lambda); arb_clear(one_minus);
  arb_clear(a); arb_clear(b); arb_clear(c); arb_clear(d);
}

// Table 1: certify the first positive zero of each q_n, for 2<=n<=21.
// Store its enclosure in roots[n], indexed by half the number of polygon sides.
void certify_polygon_roots(arb_ptr roots, slong prec){
  arb_t left, right, res, pi, bound;
  fmpq_t q;
  arb_init(left); arb_init(right); arb_init(res);
  arb_init(pi); arb_init(bound);
  fmpq_init(q);

  // Lemma 4.3 uses sqrt(pi/2)<4/3 to give positivity on (0,3/2].
  arb_const_pi(pi,prec); arb_div_ui(res,pi,2,prec); arb_sqrt(res,res,prec);
  arb_set_ui(bound,4); arb_div_ui(bound,bound,3,prec);
  check_less(res,bound,"small radius reduction");

  ulong total_boxes = 0;
  slong max_depth = 0;
  for(ulong j = 0; j < number_polygon_root_brackets; j++){
    ulong n = polygon_root_brackets[j].n;
    fmpq_set_str(q,polygon_root_brackets[j].lower_exact,10);
    fmpq_canonicalise(q);
    arb_set_fmpq(left,q,prec);
    fmpq_set_str(q,polygon_root_brackets[j].upper_exact,10);
    fmpq_canonicalise(q);
    arb_set_fmpq(right,q,prec);
    arb_polygon_q(res,left,n,prec);
    check_positive(res,"left root bracket",n);
    arb_polygon_q(res,right,n,prec);
    check_negative(res,"right root bracket",n);

    // Check q_n>0 on [3/2,lower] to rule out earlier zeros.
    subdivision_stats stats;
    if(!test_positivity(arb_polygon_q_interval,&n,"3/2",
                        polygon_root_brackets[j].lower_exact,120,prec,stats))
      failed("positivity before root",n);
    total_boxes += stats.boxes;
    max_depth = max(max_depth,stats.depth);
    arb_union(roots+n,left,right,prec);
  }

  if(print_output)
    cout << "POLYGON ROOTS: PASSED (" << total_boxes
         << " boxes, depth " << max_depth << ")" << endl;

  arb_clear(left); arb_clear(right); arb_clear(res);
  arb_clear(pi); arb_clear(bound);
  fmpq_clear(q);
}

// Proposition 4.4: check the eight inequalities in (4.7).
void certify_finite_chain(arb_srcptr roots, const arb_t j11){
  arb_srcptr chain[9] = {roots+6, roots+7, roots+8, roots+9, j11,
                         roots+5, roots+4, roots+3, roots+2};
  const char *names[9] = {"P12", "P14", "P16", "P18", "disk",
                           "P10", "P8", "P6", "P4"};
  for(int j = 0; j < 8; j++){
    if(!arb_gt(chain[j],chain[j+1])){
      cerr << names[j] << ">" << names[j+1] << ": FAILED" << endl;
      exit(1);
    }
  }
  passed("FINITE COMPARISON");
}

// Proposition 4.5: check |roots[n]-j11-C0/n^6| < 4/n^8 for 9<=n<=20.
// roots, j11 and C0 have already been certified.
void certify_bridge(arb_srcptr roots, const arb_t j11, const arb_t C0, slong prec){
  arb_t center, error, res;
  arb_init(center); arb_init(error); arb_init(res);

  for(ulong n = 9; n <= 20; n++){
    arb_div_ui(center,C0,n_pow(n,6),prec);
    arb_add(center,center,j11,prec);
    arb_set_ui(error,4); arb_div_ui(error,error,n_pow(n,8),prec);
    arb_sub(res,roots+n,center,prec); arb_add(res,res,error,prec);
    check_positive(res,"lower bridge bound",n);
    arb_sub(res,center,roots+n,prec); arb_add(res,res,error,prec);
    check_positive(res,"upper bridge bound",n);
  }

  passed("BRIDGE ESTIMATE");
  arb_clear(center); arb_clear(error); arb_clear(res);
}

// Lemmas 4.8, 4.11, Proposition 4.10 and Theorem 4.6.
// Check the bounds at n=21; the paper extends them to n>=21.
void certify_effective_bounds(const arb_t j11, const arb_t C0, slong prec){
  ulong n = 21;
  arb_t pi, c0677, a, b, c, bound, delta, R_sq, R, symmetric_difference, base;
  arb_init(pi); arb_init(c0677); arb_init(a); arb_init(b); arb_init(c);
  arb_init(bound); arb_init(delta); arb_init(R_sq); arb_init(R);
  arb_init(symmetric_difference); arb_init(base);
  arb_const_pi(pi,prec);
  arb_set_str(c0677,"0.677",prec);

  // h_n(u/n)=G0(u)/n^2+Z_n(u)/n^4, with ||Z_n||_infinity < 0.677.
  arb_div_ui(a,pi,2,prec); arb_pow_ui(a,a,4,prec); arb_div_ui(a,a,9,prec);
  check_less(a,c0677,"bound for Z_n");
  arb_pow_ui(a,pi,3,prec);
  arb_sqrt_ui(b,3,prec); arb_mul_ui(b,b,18,prec); arb_div(a,a,b,prec);
  arb_one(bound);
  check_less(a,bound,"integral of |G_0|");

  arb_pow_ui(a,pi,5,prec); arb_div_ui(a,a,720,prec); arb_sqrt(a,a,prec);
  arb_sqrt(b,pi,prec); arb_mul(b,b,c0677,prec); arb_div_ui(b,b,n*n,prec);
  arb_add(a,a,b,prec); arb_mul(a,a,a,prec); arb_mul_ui(a,a,2,prec);
  arb_one(bound);
  check_less(a,bound,"second moment");

  arb_pow_ui(a,pi,5,prec); arb_mul(a,a,c0677,prec);
  arb_mul_ui(a,a,3,prec); arb_div_ui(a,a,720,prec);
  arb_mul(b,c0677,c0677,prec); arb_mul_ui(b,b,3,prec); arb_div_ui(b,b,n*n,prec);
  arb_pow_ui(c,c0677,3,prec); arb_mul(c,c,pi,prec);
  arb_div_ui(c,c,n*n*n*n,prec);
  arb_add(a,a,b,prec); arb_add(a,a,c,prec); arb_mul_ui(a,a,2,prec);
  arb_set_ui(bound,7); arb_div_ui(bound,bound,4,prec);
  check_less(a,bound,"third moment");

  // Check n^8*I_4 < 3/10 (Lemma 4.8).
  arb_pow_ui(a,pi,9,prec); arb_div_ui(a,a,241920,prec); arb_root_ui(a,a,4,prec);
  arb_root_ui(b,pi,4,prec); arb_mul(b,b,c0677,prec); arb_div_ui(b,b,n*n,prec);
  arb_add(a,a,b,prec); arb_pow_ui(a,a,4,prec); arb_mul_ui(a,a,2,prec);
  arb_set_ui(bound,3); arb_div_ui(bound,bound,10,prec);
  check_less(a,bound,"fourth moment");

  arb_mul(a,pi,pi,prec); arb_div_ui(a,a,12,prec);
  arb_div_ui(b,c0677,n*n,prec); arb_add(a,a,b,prec); arb_div_ui(a,a,n*n,prec);
  arb_one(bound); arb_div_ui(bound,bound,500,prec);
  check_less(a,bound,"supremum of h_n");
  // Mean error coefficient: 0.9856*(7/4)+(13/24)*(3/10)=1.8873 < 1.89.
  arb_set_str(a,"0.9856",prec); arb_mul_ui(a,a,7,prec); arb_div_ui(a,a,4,prec);
  arb_set_ui(b,13); arb_div_ui(b,b,24,prec);
  arb_mul_ui(b,b,3,prec); arb_div_ui(b,b,10,prec); arb_add(a,a,b,prec);
  arb_set_str(bound,"1.89",prec);
  check_less(a,bound,"error in a0(j11)");

  // Geometry and slope estimates. It is enough to evaluate them at n=21.
  arb_div_ui(delta,pi,2*n,prec);
  arb_mul_ui(a,delta,2,prec); arb_sin(b,a,prec); arb_div(R_sq,a,b,prec);
  arb_sqrt(R,R_sq,prec);
  arb_tan(a,delta,prec); arb_mul(symmetric_difference,pi,delta,prec);
  arb_mul(symmetric_difference,symmetric_difference,a,prec);

  arb_set_str(a,"0.1",prec); arb_add(a,a,j11,prec);
  arb_mul(a,a,R_sq,prec); arb_mul(a,a,symmetric_difference,prec);
  arb_set_str(bound,"0.07",prec);
  check_less(a,bound,"slope error");
  arb_mul(a,R_sq,symmetric_difference,prec);
  arb_set_str(bound,"0.02",prec);
  check_less(a,bound,"positive slope error");

  arb_set_str(a,"0.16",prec); arb_mul(a,a,a,prec);
  arb_set_str(bound,"0.5",prec);
  check_less(a,bound,"geometric tail ratio");
  arb_sqrt_ui(a,77,prec); arb_set_str(bound,"8.77",prec);
  check_greater(a,bound,"sqrt(77)");
  arb_mul_ui(a,pi,2,prec);
  check_greater(bound,a,"8.77 > 2*pi");

  arb_mul(a,j11,R,prec); arb_set_str(bound,"3.84",prec);
  check_less(a,bound,"j11*R_21");
  arb_set_str(a,"0.196868",prec);
  arb_set_ui(b,19); arb_div_ui(b,b,10,prec); arb_div_ui(b,b,n*n,prec);
  arb_sub(a,a,b,prec);
  check_positive(a,"q_n(j11) endpoint");

  // Fourier tail: n^8*4*pi*1.92^(2n)/(2n)! < 10^-20.
  // Here j11*R_21/2 < 1.92. The ratio of consecutive bounds is checked below.
  arb_set_str(base,"1.92",prec);
  arb_pow_ui(a,base,2*n,prec); arb_mul(a,a,pi,prec); arb_mul_ui(a,a,4,prec);
  arb_fac_ui(b,2*n,prec); arb_div(a,a,b,prec);
  arb_mul_ui(a,a,n_pow(n,8),prec);
  arb_set_str(bound,"1e-20",prec);
  check_less(a,bound,"tail in q_n(j11)");

  // Add the Fourier tail to the mean error: 1.89+10^-20 < 19/10.
  arb_set_str(a,"1.89",prec); arb_add(a,a,bound,prec);
  arb_set_ui(bound,19); arb_div_ui(bound,bound,10,prec);
  check_less(a,bound,"error in q_n(j11)");

  arb_set_ui(a,n_pow(n+1,8)); arb_div_ui(a,a,n_pow(n,8),prec);
  arb_mul(a,a,base,prec); arb_mul(a,a,base,prec);
  arb_div_ui(a,a,(2*n+1)*(2*n+2),prec);
  arb_one(bound);
  check_less(a,bound,"decreasing tail in q_n(j11)");

  arb_mul(base,j11,R,prec); arb_div_ui(base,base,2,prec);
  arb_pow_ui(a,base,2*n-1,prec); arb_mul(a,a,pi,prec);
  arb_mul(a,a,R,prec); arb_mul_ui(a,a,4,prec);
  arb_fac_ui(b,2*n-1,prec); arb_div(a,a,b,prec);
  arb_one(bound); arb_div_ui(bound,bound,n_pow(n,4),prec);
  check_less(a,bound,"tail in the derivative");

  arb_set_ui(a,n_pow(n+1,4)); arb_div_ui(a,a,n_pow(n,4),prec);
  arb_mul(a,a,base,prec); arb_mul(a,a,base,prec);
  arb_div_ui(a,a,(2*n)*(2*n+1),prec);
  arb_one(bound);
  check_less(a,bound,"decreasing tail in the derivative");

  // Effective zero estimate: |q_n(j11+C0/n^6)| < 1.93/n^8.
  // Check 1.93/0.53<4 and C0/n^6+4/n^8<0.1 to stay in the slope interval.
  arb_set_ui(a,19); arb_div_ui(a,a,10,prec);
  arb_mul_ui(b,C0,2,prec); arb_div_ui(b,b,n*n,prec); arb_add(a,a,b,prec);
  arb_mul(b,C0,C0,prec); arb_mul(b,b,pi,prec); arb_mul(b,b,R_sq,prec);
  arb_div_ui(b,b,4*n*n*n*n,prec); arb_add(a,a,b,prec);
  arb_set_str(bound,"1.93",prec);
  check_less(a,bound,"effective error");
  arb_set_str(a,"1.93",prec); arb_set_str(b,"0.53",prec); arb_div(a,a,b,prec);
  arb_set_ui(bound,4);
  check_less(a,bound,"1.93/0.53");
  arb_div_ui(a,C0,n_pow(n,6),prec);
  arb_set_ui(b,4); arb_div_ui(b,b,n_pow(n,8),prec); arb_add(a,a,b,prec);
  arb_set_str(bound,"0.1",prec);
  check_less(a,bound,"interval for the effective estimate");

  // Theorem 1.5, Section 4.4: exact N_9, N_10, N_11 and the n=12 check.
  ulong expected[3] = {279410415100UL,762329564900UL,1864174692448UL};
  for(ulong k = 9; k <= 11; k++){
    ulong value = 149*k*k*(k+1)*(k+1)
                  *(n_pow(k+1,6)-n_pow(k,6))
                  -2000*(n_pow(k+1,8)+n_pow(k,8));
    if(value != expected[k-9])
      failed("integer check N",k);
  }
  if(894UL*n_pow(12,8) <= 4000UL*n_pow(13,7))
    failed("last integer check");

  passed("EFFECTIVE ESTIMATE");
  arb_clear(pi); arb_clear(c0677); arb_clear(a); arb_clear(b); arb_clear(c);
  arb_clear(bound); arb_clear(delta); arb_clear(R_sq); arb_clear(R);
  arb_clear(symmetric_difference); arb_clear(base);
}

int main(int argc, char *argv[]){
  slong prec = 384;
  for(int j = 1; j < argc; j++){
    if(strcmp(argv[j],"--quiet") == 0) print_output = 0;
    else if(strcmp(argv[j],"--bits") == 0 && j+1 < argc) prec = atol(argv[++j]);
    else if(strcmp(argv[j],"--help") == 0){
      cout << "Usage: certify [--bits N] [--quiet]" << endl;
      return 0;
    }
    else{
      cerr << "Unknown option: " << argv[j] << endl;
      return 1;
    }
  }
  if(prec < 256){
    cerr << "The precision must be at least 256 bits" << endl;
    return 1;
  }

  arb_t j11, C0;
  arb_init(j11); arb_init(C0);
  slong number_roots = (slong)polygon_root_brackets[number_polygon_root_brackets-1].n+1;
  arb_ptr roots = _arb_vec_init(number_roots);

  certify_constants(j11,C0,prec);
  certify_bessel_bounds(j11,prec);
  certify_elementary_bounds(prec);
  certify_master(prec);
  certify_polygon_roots(roots,prec);
  certify_finite_chain(roots,j11);
  certify_bridge(roots,j11,C0,prec);
  certify_effective_bounds(j11,C0,prec);

  if(print_output) cout << "CERTIFICATE PASSED" << endl;

  arb_clear(j11); arb_clear(C0);
  _arb_vec_clear(roots,number_roots);
  flint_cleanup();
  return 0;
}
