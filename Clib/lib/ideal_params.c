#include "gen_param.h"
#include <fmpz.h>
#include <fmpz_mat.h>
#include <flint.h>

/*B should be small*/
fmpz_t *
split_primes_below_B(fmpz_t f, ulong *nb_primes, ulong d0, ulong B)
{
    fmpz_t *split_primes = flint_malloc(sizeof(fmpz_t) * B);
    if (!split_primes)
        return NULL;

    fmpz_t p, tmp;
    fmpz_init_set_ui(p, 2);
    fmpz_init(tmp);

    ulong j = 0;
    while (fmpz_cmp_ui(p, B) < 0)
    {
        fmpz_nextprime(p, p, 0);
        if (fmpz_mod_ui(tmp, p, 4) != 1)
            continue;

        /*p is split in Z[i]*/
        fmpz_init_set(split_primes[j], p);
        j++;
    }
    split_primes = flint_realloc(split_primes, (j) * sizeof(fmpz_t));
    *nb_primes = j;

    return split_primes;
}

/*The primes should be congruent to 1 mod 4*/
fmpz_t *
ordered_sqrt_minus_one_mod_l(fmpz_t *primes, ulong nb_primes)
{
    fmpz_t *sqrts_of_minus_one_mod_l = flint_malloc(sizeof(fmpz_t) * nb_primes);
    if (!sqrts_of_minus_one_mod_l)
        return NULL;

    fmpz_t l, _1;
    fmpz_init(l);
    fmpz_init_set_si(_1, -1);

    for (ulong j = 0; j < nb_primes; j++)
    {
        fmpz_init(sqrts_of_minus_one_mod_l[j]);
        fmpz_set(l, primes[j]);

        fmpz_sqrtmod(sqrts_of_minus_one_mod_l[j], _1, l);
    }

    return sqrts_of_minus_one_mod_l;
}
