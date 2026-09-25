// Implementation file
// Checks strict bounds on Arb balls and functions over rational intervals.

#include "methods.h"

#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include "flint/fmpq.h"

using namespace std;

void failed(const char *name){
  cerr << name << ": FAILED" << endl;
  exit(1);
}

void failed(const char *name, ulong n){
  cerr << name << ", n=" << n << ": FAILED" << endl;
  exit(1);
}

// A comparison passes only when Arb proves it for the entire balls.
// An inconclusive enclosure terminates the certificate just like a false bound.
void check_positive(const arb_t x, const char *name){
  if(!arb_is_positive(x)) failed(name);
}

void check_positive(const arb_t x, const char *name, ulong n){
  if(!arb_is_positive(x)) failed(name,n);
}

void check_negative(const arb_t x, const char *name){
  if(!arb_is_negative(x)) failed(name);
}

void check_negative(const arb_t x, const char *name, ulong n){
  if(!arb_is_negative(x)) failed(name,n);
}

void check_less(const arb_t x, const arb_t y, const char *name){
  if(!arb_lt(x,y)) failed(name);
}

void check_less(const arb_t x, const arb_t y, const char *name, ulong n){
  if(!arb_lt(x,y)) failed(name,n);
}

void check_greater(const arb_t x, const arb_t y, const char *name){
  if(!arb_gt(x,y)) failed(name);
}

void check_greater(const arb_t x, const arb_t y, const char *name, ulong n){
  if(!arb_gt(x,y)) failed(name,n);
}

void check_inside(const arb_t x, const arb_t left, const arb_t right,
                  const char *name){
  check_greater(x,left,name);
  check_less(x,right,name);
}

enum check_type { POSITIVE, LOWER_BOUND, ABSOLUTE_BOUND };

int check_interval_recursive(arb_interval_func F, const void *params,
                             const fmpq_t left, const fmpq_t right,
                             const arb_t bound, check_type type, slong level,
                             slong max_depth, slong prec,
                             subdivision_stats& stats){
  arb_t a, b, I, res, abs_res;
  arb_init(a); arb_init(b); arb_init(I); arb_init(res); arb_init(abs_res);
  arb_set_fmpq(a,left,prec);
  arb_set_fmpq(b,right,prec);
  // I encloses [left,right]; F is evaluated on the whole interval.
  arb_union(I,a,b,prec);
  F(res,I,params,prec);

  // Accept only a strict bound; overlapping balls require subdivision.
  int success = 0;
  if(type == POSITIVE) success = arb_is_positive(res);
  if(type == LOWER_BOUND) success = arb_gt(res,bound);
  if(type == ABSOLUTE_BOUND){
    arb_abs(abs_res,res);
    success = arb_lt(abs_res,bound);
  }

  if(success){
    stats.boxes++;
    stats.depth = max(stats.depth,level);
    arb_clear(a); arb_clear(b); arb_clear(I); arb_clear(res); arb_clear(abs_res);
    return 1;
  }
  if(level == max_depth){
    cout << "Max depth achieved on ";
    arb_printd(I,20);
    cout << endl;
    arb_clear(a); arb_clear(b); arb_clear(I); arb_clear(res); arb_clear(abs_res);
    return 0;
  }

  arb_clear(a); arb_clear(b); arb_clear(I); arb_clear(res); arb_clear(abs_res);
  
  fmpq_t middle;
  fmpq_init(middle);
  fmpq_add(middle,left,right);
  fmpq_div_2exp(middle,middle,1);
  int result_left = check_interval_recursive(F,params,left,middle,bound,type,
                                              level+1,max_depth,prec,stats);
  int result_right = 0;
  if(result_left)
    result_right = check_interval_recursive(F,params,middle,right,bound,type,
                                            level+1,max_depth,prec,stats);
  fmpq_clear(middle);
  return result_left && result_right;
}

int test_interval(arb_interval_func F, const void *params,
                  const char *left_string, const char *right_string,
                  const arb_t bound, check_type type, slong max_depth,
                  slong prec, subdivision_stats& stats){
  fmpq_t left, right;
  fmpq_init(left); fmpq_init(right);
  fmpq_set_str(left,left_string,10);
  fmpq_canonicalise(left);
  fmpq_set_str(right,right_string,10);
  fmpq_canonicalise(right);
  if(fmpq_cmp(left,right) >= 0){
    cerr << "Empty interval in subdivision" << endl;
    exit(1);
  }
  stats.boxes = 0; stats.depth = 0;
  int result = check_interval_recursive(F,params,left,right,bound,type,0,
                                        max_depth,prec,stats);
  fmpq_clear(left); fmpq_clear(right);
  return result;
}

int test_positivity(arb_interval_func F, const void *params,
                    const char *left, const char *right,
                    slong max_depth, slong prec, subdivision_stats& stats){
  arb_t dummy;
  arb_init(dummy); arb_zero(dummy);
  int result = test_interval(F,params,left,right,dummy,POSITIVE,
                             max_depth,prec,stats);
  arb_clear(dummy);
  return result;
}

int test_lower_bound(arb_interval_func F, const void *params,
                     const char *left, const char *right, const arb_t bound,
                     slong max_depth, slong prec, subdivision_stats& stats){
  return test_interval(F,params,left,right,bound,LOWER_BOUND,
                       max_depth,prec,stats);
}

int test_absolute_bound(arb_interval_func F, const void *params,
                        const char *left, const char *right, const arb_t bound,
                        slong max_depth, slong prec, subdivision_stats& stats){
  return test_interval(F,params,left,right,bound,ABSOLUTE_BOUND,
                       max_depth,prec,stats);
}
