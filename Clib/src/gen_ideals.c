#include "error.h"
#include "gen_param.h"

#include <fmpz.h>
#include <fmpz_mat.h>
#include <flint.h>

static void
print_data(fmpz_t *split_primes, fmpz_t *sqrts_mod_l, ulong nb_primes, FILE* F);

static void
fprint_and_clear_data(fmpz_t *split_primes, fmpz_t *sqrts_mod_l, ulong nb_primes, FILE* F);


int main(int argc, char *argv[]){
    if(argc != 3){
        flint_printf("Usage:\n\t -Input: <B> <fichier> \n\t -Output: liste de premiers l_i décomposés dans Z[i] et racines carrées de -1 mod l_i\n\n");
        exit(1);
    }

    FILE *F = fopen(argv[2], "w");
    if(!F) exit(1);

    fmpz_t B;
    fmpz_init(B);
    
    if(fmpz_set_str(B, argv[1], 10) == -1){
        fclose(F);
        fmpz_clear(B);
        exit(1);
    }
    if(fmpz_cmp_ui(B, 10000)>0){
        fclose(F);
        fmpz_clear(B);
        flint_printf("Slow down, try another B\n");
        exit(1);
    }
    fmpz_t f;
    fmpz_init(f);
    
    /*B is assumed to be 2**FLINT_BITS-1 long at most*/
    ulong nb_primes = 0;
    fmpz_t *split_primes = split_primes_below_B(f, &nb_primes, -4, fmpz_get_ui(B));
    fmpz_clear(B);
    fmpz_clear(f);
    fmpz_t *sqrts_mod_l = ordered_sqrt_minus_one_mod_l(split_primes, nb_primes);


    print_data(split_primes, sqrts_mod_l, nb_primes, F);
    fprint_and_clear_data(split_primes, sqrts_mod_l, nb_primes, F);

    return 0;
}

static void
print_data(fmpz_t *split_primes, fmpz_t *sqrts_mod_l, ulong nb_primes, FILE* F)
{
    for(ulong i = 0; i<nb_primes-1; i++){
        flint_printf("(");
        fmpz_print(sqrts_mod_l[i]);
        flint_printf(")**2=-1 mod (");
        fmpz_print(split_primes[i]);
        flint_printf(")\n");
    }
    flint_printf("(");

    fmpz_print(sqrts_mod_l[nb_primes-1]);
    flint_printf(")**2=-1 mod (");
    fmpz_print(split_primes[nb_primes-1]);
    flint_printf(")\n");
}

static void
fprint_and_clear_data(fmpz_t *split_primes, fmpz_t *sqrts_mod_l, ulong nb_primes, FILE* F)
{
    fmpz_mat_t M;
    fmpz_mat_init(M, 2, 2);

    fmpz_t n;
    fmpz_init_set_ui(n, nb_primes);
    fmpz_fprint(F, n);
    fputc('\n',F);
    fmpz_clear(n);

    fmpz arr[4];
    for(int i = 0; i<4; i++) fmpz_init_set_ui(arr+i, 0);

    for(ulong i = 0; i<nb_primes; i++){
        fmpz_set(arr, split_primes[i]);
        //fmpz_sub(arr+1, split_primes[i], sqrts_mod_l[i]);
        fmpz_neg(arr+1, sqrts_mod_l[i]);
        fmpz_set_ui(arr+3, 1);

        fmpz_clear(split_primes[i]);
        fmpz_clear(sqrts_mod_l[i]);

        fmpz_arr_to_mat(M, arr, 4, 2, 2);
        fmpz_mat_fprint(F, M);
        fputc('\n', F);
    }

    for(int i = 0; i<4; i++) fmpz_clear(arr+i);

    fmpz_mat_clear(M);
    flint_free(split_primes);
    flint_free(sqrts_mod_l);
    fclose(F);
}