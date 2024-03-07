#include "gen_param.h"
#include <fmpz.h>
#include <fmpz_vec.h>
#include <fmpz_mat.h>


int main(int argc, char *argv[])
{
    if(argc!=2){
        flint_printf("Usage: ./find_conductor <filename>\n");
        exit(1);
    }

    FILE *f = fopen(argv[1], "r");
    if(!f){
        flint_printf("Invalid file\n");
        exit(1);
    }

    fmpz_t nb_primes;
    fmpz_init(nb_primes);
    if(!fmpz_fread(f, nb_primes))
    {
        fmpz_clear(nb_primes);
        fclose(f);
        flint_printf("Invalid number format in file\n");
        exit(1);
    }

    fmpz_mat_t I, J, M;
    fmpz_mat_init(I, 2, 2);
    fmpz_mat_init(J, 2, 2);
    fmpz_mat_init(M, 4, 2);

    if(!fmpz_mat_fread(f, J))
    {
        flint_printf("Invalid matrix format in file\n");
        fmpz_mat_clear(I);
        fmpz_mat_clear(J);
        fmpz_mat_clear(M);
        fmpz_clear(nb_primes);
        exit(1);
    }
    for(ulong i = 0; i < fmpz_get_ui(nb_primes)-1; i++)
    {
        if(!fmpz_mat_fread(f, I))
        {
            flint_printf("Invalid matrix format in file\n");
            fmpz_mat_clear(I);
            fmpz_mat_clear(J);
            fmpz_mat_clear(M);
            fmpz_clear(nb_primes);
            exit(1);
        }

        // if(!fmpz_mat_fread(f, J))
        // {
        //     flint_printf("Invalid matrix format in file\n");
        //     fmpz_mat_clear(I);
        //     fmpz_mat_clear(J);
        //     fmpz_mat_clear(M);
        //     fmpz_clear(nb_primes);
        //     exit(1);
        // }

        ideal_product(J, M, I, J);
        flint_printf("\nIJ=");
        fmpz_mat_print_pretty(J);
        flint_printf("\n");
    }
    

    fmpz_mat_clear(I);
    fmpz_mat_clear(J);
    fmpz_mat_clear(M);
    fmpz_clear(nb_primes);
    exit(0);
}