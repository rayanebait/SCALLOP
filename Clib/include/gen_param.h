#ifndef GEN_PARAM_H
#define GEN_PARAM_H

#include <fmpz.h>
#include <fmpz_mat.h>

fmpz_t *
split_primes_below_B(fmpz_t f, ulong *nb_primes, ulong d0, ulong B);

fmpz_t * 
ordered_sqrt_minus_one_mod_l(fmpz_t* primes, ulong nb_primes);


int fmpz_arr_to_mat(fmpz_mat_t M, fmpz *arr, ulong len_arr, ulong n, ulong m);

void ideal_to_hnf(fmpz_mat_t I, fmpz *entries);

void ideal_product(fmpz_mat_t IJ, fmpz_mat_t M, fmpz_mat_t I, fmpz_mat_t J);

#endif 