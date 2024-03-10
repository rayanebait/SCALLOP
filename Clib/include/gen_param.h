#ifndef GEN_PARAM_H
#define GEN_PARAM_H

#include <fmpz.h>
#include <fmpz_mat.h>
#include "qfb.h"

fmpz_t *
split_primes_below_B(fmpz_t f, ulong *nb_primes, ulong d0, ulong B);

fmpz_t *
ordered_sqrt_minus_one_mod_l(fmpz_t *primes, ulong nb_primes);

int fmpz_arr_to_mat(fmpz_mat_t M, fmpz *arr, ulong len_arr, ulong n, ulong m);

void ideal_to_hnf(fmpz_mat_t I, fmpz *entries);

void ideal_product(fmpz_mat_t IJ, fmpz_mat_t M, fmpz_mat_t I, fmpz_mat_t J);

void ideal_to_qfb(qfb_t q, fmpz_t r, fmpz_t f);

void prime_ideals_to_qfbs(qfb *forms, ulong nb_primes, fmpz_t f, fmpz *r);

void ideal_to_qfb(qfb_t q, fmpz_t r, fmpz_t f);

void qfb_a_bound(fmpz_t L, fmpz_t D);

void qfb_discrete_log_prime_order(
    fmpz_t pdlog, qfb_t h, qfb_t g,
    fmpz_t p, fmpz_t n_1, fmpz_t n_2,
    fmpz_t D, fmpz_t L);

void qfb_discrete_log_prime_power(
    fmpz_t ppdlog, qfb_t h, qfb_t g,
    fmpz_t prime_power, fmpz_t prime,
    ulong power, fmpz_t n_1, fmpz_t n_2,
    fmpz_t D, fmpz_t L);

void qfb_discrete_log(fmpz_t dlog, qfb_t h, qfb_t g, fmpz_t f, fmpz *factors, ulong nb_factors);

#endif