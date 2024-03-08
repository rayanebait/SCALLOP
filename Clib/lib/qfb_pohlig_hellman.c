#include <fmpz.h>
#include <fmpz_mat.h>
#include "qfb.h"

/*Only for I prime in Z[i] with I=[1+r^2, f(i-r)] and r^2=-1 mod N(I)*/
void ideal_to_qfb(qfb_t q, fmpz_t r, fmpz_t f)
{
    fmpz_t a, b, c;
    fmpz_init_ui(a, 1);
    fmpz_init(b, r);
    fmpz_init(c, f);

    fmpz_addmul(a, b, b);
    fmpz_mul(b, b, c);
    fmpz_mul_ui(b, b, 2);

    fmpz_mul(c, f, f);

    qfb_init(q);
    qfb_setcoeffs(q, a, b, c);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}

/*Pollard's original method where alpha is a
 generator, and beta is anything. And with a
 partition on the reduced forms given by:
        -S_0={1<=a<=n_1},
        -S_1={n_1<a<= n_2},
        -S_2={n_2<a<=sqrt(D/3)}
*/
static void
next_qfb(qfb_t g, fmpz_t a, fmpz_t b,
         qfb_t alpha, qfb_t beta, fmpz_t n_1,
         fmpz_t n_2, fmpz_t D, fmpz_t L)
{
    qfb_reduce(g, g, D);

    if (fmpz_cmp(g->a, n_1) <= 0)
    {
        qfb_nucomp(g, g, alpha, D, L);
        fmpz_add_ui(a, a, 1);
    }
    else if (fmpz_cmp(g->a, n_2) <= 0)
    {
        qfb_nudupl(g, g, D, L);
        fmpz_mul_ui(a, a, 2);
        fmpz_mul_ui(b, b, 2);
    }
    else
    {
        qfb_nucomp(g, g, beta, D, L);
        fmpz_add_ui(a, a, 1);
    }
}

/*Currently uses rho-Pollard*/
qfb_discrete_log_prime_order(fmpz_t pdlog, qfb_t h, qfb_t g, fmpz_t p)
{
    qfb_t y, y2;
    qfb_init(y);
    qfb_init(y2);

    qfb_set(y1, h);
    qfb_set(y2, g);

    while (qfb_equal(y, y2) != 1)
    {
    }
}

void prime_ideals_to_qfbs(qfb *forms, fmpz_t f, fmpz *r)
{
}
/*Computes Pohlig hellman reduction, factors should be ordered*/
void qfb_discrete_log(fmpz_t dlog, qfb_t h, qfb_t g, fmpz_t f, fmpz *factors, ulong nb_factors)
{
    fmpz_t g_i;
    fmpz_init_set(g_i, factors);
    for (ulong i = 1; i < nb_factors; i++)
    {
    }
}
