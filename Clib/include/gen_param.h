#ifndef GEN_PARAM_H
#define GEN_PARAM_H

#include <fmpz.h>
#include <fmpz_mat.h>
#include "qfb.h"

fmpz_t *
split_primes_below_B(fmpz_t f, ulong *nb_primes, ulong d0, ulong B);

fmpz_t * 
ordered_sqrt_minus_one_mod_l(fmpz_t* primes, ulong nb_primes);


int fmpz_arr_to_mat(fmpz_mat_t M, fmpz *arr, ulong len_arr, ulong n, ulong m);

void ideal_to_hnf(fmpz_mat_t I, fmpz *entries);

void ideal_product(fmpz_mat_t IJ, fmpz_mat_t M, fmpz_mat_t I, fmpz_mat_t J);


void qfb_init_set(qfb_t q, fmpz_t a, fmpz_t b, fmpz_t c);

void ideal_to_qfb(qfb_t q, fmpz_t r, fmpz_t f);
void qfb_discrete_log_prime_order(fmpz_t dlog_p, qfb_t h, qfb_t g, fmpz_t p);

void qfb_discrete_log(fmpz_t dlog, qfb_t h, qfb_t g, fmpz_t f, fmpz *factors, ulong nb_factors);

#endif 