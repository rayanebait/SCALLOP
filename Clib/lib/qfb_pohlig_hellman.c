#include <fmpz.h>
#include <fmpz_mat.h>
#include "qfb.h"
#include "gen_param.h"
#include <unistd.h>

/*Only for I prime in Z[i] with I=[1+r^2, f(i-r)] and r^2=-1 mod N(I)*/
void ideal_to_qfb(qfb_t q, fmpz_t r, fmpz_t f)
{
    fmpz_t a, b, c;
    fmpz_init_set_ui(a, 1);
    fmpz_init_set(b, r);
    fmpz_init_set(c, f);

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

void
prime_ideals_to_qfbs(qfb *forms, ulong nb_primes, fmpz_t f, fmpz *r)
{
    fmpz_t D;
    fmpz_init(D);
    fmpz_mul(D, f,f);
    fmpz_mul_si(D, D,-4);

    for(ulong i = 0; i<nb_primes; i++){
        ideal_to_qfb(forms+i, r+i, f);
        qfb_reduce(forms+i, forms+i, D);
    }
    fmpz_clear(D);
}

/*Pollard's original method where g is a
 generator, and h is anything. And with a
 partition on the reduced forms given by:
        -S_0={1<=a<=n_1},
        -S_1={n_1<a<= n_2},
        -S_2={n_2<a<=sqrt(D/3)}
    for n_1=L/3, n_2=2*n_1, L=sqrt(D/3)
*/

/*y=g^a*h^b*/
void next_qfb(qfb_t y, fmpz_t a, fmpz_t b,
              qfb_t g, qfb_t h, fmpz_t n_1,
              fmpz_t n_2, fmpz_t D, fmpz_t L)
{
    qfb_reduce(y, y, D);
    if (fmpz_cmp(y->a, n_1) <= 0)
    {
        qfb_nucomp(y, y, g, D, L);
        fmpz_add_ui(a, a, 1);
    }
    else if (fmpz_cmp(y->a, n_2) <= 0)
    {
        // flint_printf("case 2\n");
        qfb_nudupl(y, y, D, L);
        fmpz_mul_ui(a, a, 2);
        fmpz_mul_ui(b, b, 2);
    }
    else
    {
        qfb_nucomp(y, y, h, D, L);
        fmpz_add_ui(b, b, 1);
    }

    qfb_reduce(y,y,D);
}

void next_2_qfb(qfb_t y, fmpz_t a, fmpz_t b,
                qfb_t g, qfb_t h, fmpz_t n_1,
                fmpz_t n_2, fmpz_t D, fmpz_t L)
{
    next_qfb(y, a, b, g, h, n_1, n_2, D, L);
    next_qfb(y, a, b, g, h, n_1, n_2, D, L);
}
/*Currently uses rho-Pollard*/
void qfb_discrete_log_prime_order(
    fmpz_t pdlog, qfb_t h, qfb_t g,
    fmpz_t p, fmpz_t n_1, fmpz_t n_2,
    fmpz_t D, fmpz_t L)
{
    if(qfb_equal(h,g)==1){
        fmpz_set_ui(pdlog, 1);
        return;
    }

    /*can reuse pdlog*/
    fmpz_t a1, b1, a2, b2;
    fmpz_init_set_ui(a1, 1);
    fmpz_init_set_ui(b1, 1);

    fmpz_init_set_ui(a2, 1);
    fmpz_init_set_ui(b2, 1);

    qfb_t y, y2;
    qfb_init(y);
    qfb_init(y2);

    qfb_nucomp(y, g, h, D, L);
    qfb_set(y2, y);

    next_qfb(y, a1, b1, g, h, n_1, n_2, D, L);
    next_2_qfb(y2, a2, b2, g, h, n_1, n_2, D, L);

    while (qfb_equal(y, y2) != 1)
    {
        next_qfb(y, a1, b1, g, h, n_1, n_2, D, L);
        next_2_qfb(y2, a2, b2, g, h, n_1, n_2, D, L);
        fmpz_mod(a1, a1, p);
        fmpz_mod(a2, a2, p);
        fmpz_mod(b1, b1, p);
        fmpz_mod(b2, b2, p);
    }

    /*should maybe verify for equality of b1, b2*/
    fmpz_sub(a1, a1, a2);
    fmpz_sub(b1, b2, b1);
    fmpz_invmod(b1, b1, p);
    fmpz_mul(pdlog, a1, b1);
    fmpz_mod(pdlog, pdlog, p);

    qfb_clear(y);
    qfb_clear(y2);
    fmpz_clear(a1);
    fmpz_clear(b1);
    fmpz_clear(b2);
    fmpz_clear(a2);
}

void qfb_discrete_log_prime_power(
    fmpz_t ppdlog, qfb_t h, qfb_t g,
    fmpz_t prime_power, fmpz_t prime,
    ulong power, fmpz_t n_1, fmpz_t n_2,
    fmpz_t D, fmpz_t L)
{
    if(power==1){
        qfb_discrete_log_prime_order(ppdlog, h, g, prime, n_1, n_2, D, L);
        return;
    }

    fmpz_zero(ppdlog);

    fmpz_t p_i, x_i;
    fmpz_init(x_i);
    fmpz_init_set_ui(p_i, 1);

    qfb_t g_i, h_i, g_alpha_i;
    qfb_init(g_i);
    qfb_init(h_i);
    qfb_init(g_alpha_i);
    qfb_set(g_alpha_i, g);

    for (ulong i = 1; i <= power; i++)
    {
        qfb_pow_ui(g_i, g_alpha_i, D, power - i);
        qfb_pow_ui(h_i, h, D, power - i);

        qfb_discrete_log_prime_order(x_i, h_i, g_i, prime, n_1, n_2, D, L);

        fmpz_mul(x_i, x_i, p_i);

        qfb_pow(g_i, g, D, x_i);
        qfb_nucomp(g_alpha_i, g_alpha_i, g_i, D, L);

        fmpz_add(ppdlog, ppdlog, x_i);
        fmpz_mul(p_i, p_i, prime);
    }
    fmpz_clear(x_i);
    fmpz_clear(p_i);

    qfb_clear(g_i);
    qfb_clear(h_i);
    qfb_clear(g_alpha_i);
}


void
qfb_a_bound(fmpz_t L, fmpz_t D)
{
    fmpz_fdiv_q_ui(L, D, 3);
    fmpz_abs(L, L);
    fmpz_sqrt(L, L);
}

void fmpz_pp_prod(fmpz_t pp_prod, fmpz *prime_powers, ulong nb_powers)
{
    fmpz_one(pp_prod);
    for (ulong i = 0; i < nb_powers; i++)
    {
        fmpz_mul(pp_prod, pp_prod, prime_powers + i);
    }
}

/*assumes indexing is correct and pp_prod is the product of
the prime powers without the j'th prime power*/
void fmpz_n_pp_prod(fmpz_t pp_prod, fmpz *prime_powers, ulong i, ulong j)
{
    if (i == j)
    {
        fmpz_divexact(pp_prod, pp_prod, prime_powers + i);
        return;
    }

    fmpz_divexact(pp_prod, pp_prod, prime_powers + i);
    fmpz_mul(pp_prod, pp_prod, prime_powers + j);
}

/*Computes Pohlig hellman reduction, factors should be ordered*/
void qfb_discrete_log(fmpz_t dlog, qfb_t h, qfb_t g, fmpz_t f, fmpz *factors, ulong nb_factors)
{
    fmpz_t p_i;
    fmpz_init(p_i);
    fmpz prime_powers[nb_factors];
    fmpz primes[nb_factors];
    ulong powers[nb_factors];

    /*Number of prime powers*/
    ulong j = 0;
    /*Compute the prime powers*/
    for (ulong i = 0; i < nb_factors; i++)
    {
        fmpz_set(p_i, factors + i);
        fmpz_init(prime_powers + j);
        fmpz_init_set(primes + j, p_i);
        powers[j] = 1;

        while(i < nb_factors-1 && fmpz_cmp(p_i, factors + i+1) == 0)
        {
            fmpz_set(p_i,factors+i+1);
            powers[j]++;
            i++;
        }
        fmpz_pow_ui(prime_powers + j, primes + j, powers[j]);
        j++;
    }

    fmpz_t D, L;
    fmpz_init(D);
    fmpz_init(L);
    qfb_discriminant(D, g);
    qfb_a_bound(L, D);

    fmpz_t n_1, n_2;
    fmpz_init(n_1);
    fmpz_init(n_2);

    /*naive S_0, S_1, S_2 for rho pollard, seems to work*/
    fmpz_fdiv_q_ui(n_1, L, 3);
    fmpz_mul_ui(n_2, n_1, 2);

    fmpz_t pp_prod;
    fmpz_init(pp_prod);
    fmpz_pp_prod(pp_prod, prime_powers, j);
    fmpz_n_pp_prod(pp_prod, prime_powers, 0, 0);

    qfb_t g_i, h_i;
    qfb_init(g_i);
    qfb_init(h_i);

    fmpz dlogs[j];
    for (ulong i = 0; i < j; i++)
    {
        qfb_pow(g_i, g, D, pp_prod);
        qfb_pow(h_i, h, D, pp_prod);

        fmpz_init(dlogs + i);
        qfb_discrete_log_prime_power(
            dlogs + i, h_i, g_i, prime_powers + i,
            primes + i, powers[i], n_1, n_2, D, L);

        /*set pp_prod to the product of the prime powers without i+1,
        assumes pp_prod is the product of the prime powers without
        the i'th power*/
        if(i==j-1) break;
        fmpz_n_pp_prod(pp_prod, prime_powers, i+1, i);
    }


    flint_printf("\n");
    for(int i =0; i<j;i++){
        fmpz_print(dlogs+i);
        flint_printf(" mod ");
        fmpz_print(prime_powers+i);
        flint_printf("\n");
    }
    if (fmpz_multi_CRT(dlog, prime_powers, dlogs, j, 0) == 0)
    {
        flint_printf("Something bad happened\n");
        goto cleanup;
    }
    fmpz_print(dlog);
    flint_printf("\n");

cleanup:
    qfb_clear(g_i);
    qfb_clear(h_i);

    fmpz_clear(D);
    fmpz_clear(L);
    fmpz_clear(n_1);
    fmpz_clear(n_2);
    fmpz_clear(pp_prod);

    for (ulong i = 0; i < j; i++)
    {
        fmpz_clear(prime_powers + i);
        fmpz_clear(primes + i);
        fmpz_clear(dlogs + i);
    }
    return;
}
