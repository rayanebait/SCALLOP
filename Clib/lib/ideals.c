#include "gen_param.h"

#include <fmpz.h>
#include <fmpz_mat.h>

int fmpz_arr_to_mat(fmpz_mat_t M, fmpz *arr, ulong len_arr, ulong n, ulong m)
{
    if(len_arr < n * m)
        return -1;

    
    for(ulong i = 0; i < n; i++){
        ulong ind = m * i;
        for(ulong j = 0; j < m; j++)
            fmpz_set(fmpz_mat_entry(M, i, j), arr+ind+j);
    }


    return 1;
}
/*Assumes entries has 4 elements.
 Currently inits I for convenience */
void ideal_to_hnf(fmpz_mat_t I, fmpz *entries)
{
    fmpz_mat_init(I, 2, 2);
    fmpz_arr_to_mat(I, entries, 4, 2, 2);
    fmpz_mat_hnf(I, I);
}

/*
    Computes HNF reduction of the product basis of I, J.
    Assumes I,J are 2x2 and hnf.
*/
/*M is just a 4x2 buffer, hnf reduction is performed on rows in flint*/

/*
    There is for sure much better ways to do this
*/

void ideal_product(fmpz_mat_t IJ, fmpz_mat_t M, fmpz_mat_t I, fmpz_mat_t J)
{
    fmpz_mat_t I_00_J, J_00_I;

    fmpz_mat_init(I_00_J, 2, 2);
    fmpz_mat_init(J_00_I, 2, 2);

    fmpz_mat_transpose(J_00_I, I);
    fmpz_mat_transpose(I_00_J, J);

    fmpz_mat_scalar_mul_fmpz(J_00_I, J_00_I, fmpz_mat_entry(J, 0, 0));
    fmpz_mat_scalar_mul_fmpz(I_00_J, I_00_J, fmpz_mat_entry(I, 0, 0));

    fmpz_fmms(fmpz_mat_entry(J_00_I, 0, 0),
              fmpz_mat_entry(I, 0, 1),
              fmpz_mat_entry(J, 0, 1), 
              fmpz_mat_entry(I, 1, 1),
              fmpz_mat_entry(J, 1, 1));

    fmpz_fmma(fmpz_mat_entry(J_00_I, 0, 1),
              fmpz_mat_entry(I, 0, 1),
              fmpz_mat_entry(J, 1, 1), 
              fmpz_mat_entry(I, 1, 1),
              fmpz_mat_entry(J, 0, 1));

    // fmpz_mat_print_pretty(J_00_I);
    // fmpz_mat_print_pretty(I_00_J);
    // flint_printf("\n\n");
    fmpz_mat_concat_vertical(M, I_00_J, J_00_I);
    fmpz_mat_clear(I_00_J);
    fmpz_mat_clear(J_00_I);

    // fmpz_mat_print_pretty(M);
    // flint_printf("\n\n");
    fmpz_mat_hnf(M, M);

    fmpz_mat_t window;
    fmpz_mat_window_init(window, M, 0, 0, 2, 2);
    fmpz_set(fmpz_mat_entry(IJ, 0, 0), fmpz_mat_entry(window, 1, 1));
    fmpz_set(fmpz_mat_entry(IJ, 0, 1), fmpz_mat_entry(window, 0, 1));
    fmpz_set_ui(fmpz_mat_entry(IJ, 1, 0), 0);
    fmpz_set(fmpz_mat_entry(IJ, 1, 1), fmpz_mat_entry(window, 0, 0));

    fmpz_mat_window_clear(window);
}