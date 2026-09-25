// Header file
// Methods to check signs and bounds on Arb balls and functions.

#ifndef METHODS_H_INCLUDED
#define METHODS_H_INCLUDED

#include "flint/arb.h"

void failed(const char *name);
void failed(const char *name, ulong n);

// Strict comparisons; false or inconclusive bounds terminate the certificate.
void check_positive(const arb_t x, const char *name);
void check_positive(const arb_t x, const char *name, ulong n);
void check_negative(const arb_t x, const char *name);
void check_negative(const arb_t x, const char *name, ulong n);
void check_less(const arb_t x, const arb_t y, const char *name);
void check_less(const arb_t x, const arb_t y, const char *name, ulong n);
void check_greater(const arb_t x, const arb_t y, const char *name);
void check_greater(const arb_t x, const arb_t y, const char *name, ulong n);
void check_inside(const arb_t x, const arb_t left, const arb_t right,
                  const char *name);

typedef void (*arb_interval_func)(arb_t, const arb_t, const void *, slong);

struct subdivision_stats {
  ulong boxes;
  slong depth;
};

// Return 1 only when the strict bound holds on the whole closed interval;
// return 0 if subdivision reaches max_depth without a proof.
// stats counts accepted boxes and their greatest depth (the initial box is 0).
// test_positivity checks F(x)>0.
int test_positivity(arb_interval_func F, const void *params,
                    const char *left, const char *right,
                    slong max_depth, slong prec, subdivision_stats& stats);

// Check F(x)>bound throughout [left,right].
int test_lower_bound(arb_interval_func F, const void *params,
                     const char *left, const char *right, const arb_t bound,
                     slong max_depth, slong prec, subdivision_stats& stats);

// Check |F(x)|<bound throughout [left,right].
int test_absolute_bound(arb_interval_func F, const void *params,
                        const char *left, const char *right, const arb_t bound,
                        slong max_depth, slong prec, subdivision_stats& stats);

#endif
