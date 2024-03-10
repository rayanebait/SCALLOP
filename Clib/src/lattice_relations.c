#include "error.h"
#include "gen_param.h"
#include "qfb.h"

#include <fmpz.h>
#include <fmpz_mat.h>
#include <unistd.h>

static void 
fprint_data(FILE* F, fmpz *dlogs, ulong nb_dlogs, ulong gen_ind);

static void
earlycleanup1(FILE *G, FILE *H, FILE *F,
              fmpz_t f, fmpz_t n_primes,
              fmpz_t n_factors);
static void
earlycleanup2(FILE *G, FILE *H, FILE *F,
              fmpz *sqrts, fmpz *factors,
              fmpz_t f, fmpz_t n_primes,
              fmpz_t n_factors);

static void
qfb_test_order2(fmpz_t D);

static void
single_discrete_log(
    fmpz_t dlog, qfb_t test, qfb_t g,
    qfb_t h, fmpz_t f, fmpz_t D,
    fmpz *factors, ulong nb_factors);

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        flint_printf("Usage:\n\t -Input: <fichier du conducteur> <fichier des racines> <fichier ou ecrire>\n\t -Output: Dans Cl(Z+fZ[i]), log discret de n premiers au dessus des p=1 mod 4, en base l|(5).\n\n");
        exit(1);
    }

    FILE *H = fopen(argv[1], "r");
    if (!H)
        exit(1);
    FILE *G = fopen(argv[2], "r");
    if (!G)
        exit(1);
    FILE *F = fopen(argv[3], "w");
    if (!F)
        exit(1);

    fmpz_t f, n_primes, n_factors;
    fmpz_init(n_primes);
    fmpz_init(n_factors);
    fmpz_init(f);

    if (!fmpz_fread(G, n_primes) || !fmpz_fread(H, n_factors))
    {
        flint_printf("Invalid integer format in file\n");
        earlycleanup1(G, H, F, f, n_primes, n_factors);
        exit(1);
    }

    ulong nb_primes = fmpz_get_ui(n_primes);
    ulong nb_factors = fmpz_get_ui(n_factors);

    fmpz factors[nb_factors];
    fmpz sqrts[nb_primes];
    for (ulong i = 0; i < nb_primes; i++)
    {
        fmpz_init(sqrts + i);
    }
    for (ulong i = 0; i < nb_factors; i++)
    {
        fmpz_init(factors + i);
    }

    for (ulong i = 0; i < nb_primes; i++)
    {
        if (!fmpz_fread(G, sqrts + i))
        {
            flint_printf("Invalid integer format in file\n");
            earlycleanup2(G, H, F, sqrts, factors, f, n_primes, n_factors);
            exit(1);
        }
    }
    if (!fmpz_fread(H, f))
    {
        flint_printf("Invalid integer format in file\n");
        earlycleanup2(G, H, F, sqrts, factors, f, n_primes, n_factors);
        exit(1);
    }

    for (ulong i = 0; i < nb_factors; i++)
    {
        if (!fmpz_fread(H, factors + i))
        {
            flint_printf("Invalid integer format in file\n");
            earlycleanup2(G, H, F, sqrts, factors, f, n_primes, n_factors);
            exit(1);
        }
    }

    fclose(H);
    fclose(G);
    fmpz_clear(n_primes);
    fmpz_clear(n_factors);

    qfb qfbs[nb_primes];
    prime_ideals_to_qfbs(qfbs, nb_primes, f, sqrts);
    for (ulong i = 0; i < nb_primes; i++)
    {
        fmpz_clear(sqrts + i);
    }
    fmpz_t D;
    fmpz_init(D);

    fmpz_mul(D, f, f);
    fmpz_mul_si(D, D, -4);

    qfb_t l_0;
    qfb_init(l_0);

    /*set l_0 to a prime above 5 (Should generate the class group)*/
    /*Here l_0 is of order (f-(-1/f))/2 for the 40/80 bits parameters
    so generates the class group*/

    ulong gen_ind = 0;
    qfb_set(l_0, qfbs+gen_ind);
    qfb_reduce(l_0, l_0, D);

    fmpz dlogs[nb_primes];
    for (ulong i = 0; i < nb_primes; i++)
    {
        if(i==gen_ind) continue;

        fmpz_init(dlogs + i);
        qfb_discrete_log(
            dlogs+i, qfbs+i,
            l_0, f, factors,
            nb_factors);
    }
    fprint_data(F, dlogs, nb_primes, gen_ind);

    // single_discrete_log(dlogs, qfbs, l_0, qfbs+4,f, D, factors, nb_factors);
    goto cleanup;

cleanup:
    for (ulong i = 0; i < nb_primes; i++)
    {
        qfb_clear(qfbs + i);
        fmpz_clear(sqrts + i);
        fmpz_clear(dlogs + i);
    }
    for (ulong i = 0; i < nb_factors; i++)
    {
        fmpz_clear(factors + i);
    }
    fmpz_clear(D);
    fmpz_clear(f);
    qfb_clear(l_0);
    fclose(F);
    exit(0);
}

static void 
fprint_data(FILE* F, fmpz *dlogs, ulong nb_dlogs, ulong gen_ind)
{
    flint_fprintf(F, "%wu\n", gen_ind);
    flint_fprintf(F, "%wu\n", nb_dlogs);
    for(ulong i=0; i<nb_dlogs; i++)
    {
        if(i==gen_ind) continue;
        fmpz_fprint(F, dlogs+i);
        fputc('\n', F);
    }
}

static void
single_discrete_log(
    fmpz_t dlog, qfb_t test, qfb_t g,
    qfb_t h, fmpz_t f, fmpz_t D,
    fmpz *factors, ulong nb_factors)
{

    qfb_discrete_log(
        dlog, h,
        g, f, factors,
        nb_factors);

    flint_printf("dlog of h: ");
    fmpz_print(dlog);
    flint_printf("\n");
    sleep(1);

    qfb_pow(test, g, D, dlog);
    flint_printf("g**");
    flint_printf("=\n");
    qfb_print(test);
    flint_printf("\n");
    flint_printf("h=");
    qfb_print(h);
    flint_printf("\n");
}

static void qfb_test_order2(fmpz_t D)
{
    fmpz_t a, b, c, L;

    fmpz_init_set_ui(a, 2);
    fmpz_init_set_ui(b, 2);
    fmpz_init(c);
    fmpz_set_str(c, "4544174201648970478061435463868508344033512034981", 10);
    qfb_t order2qfb;
    qfb_init(order2qfb);
    qfb_setcoeffs(order2qfb, a, b, c);

    flint_printf("Order 2 bqf:");
    qfb_print(order2qfb);
    flint_printf("\n");

    fmpz_init(L);
    qfb_a_bound(L, D);
    qfb_nudupl(order2qfb, order2qfb, D, L);

    qfb_reduce(order2qfb, order2qfb, D);

    qfb_print(order2qfb);
    qfb_clear(order2qfb);
    fmpz_clear(L);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}

static void
earlycleanup1(FILE *G, FILE *H, FILE *F, fmpz_t f, fmpz_t n_primes, fmpz_t n_factors)
{
    fclose(H);
    fclose(F);
    fclose(G);
    fmpz_clear(n_primes);
    fmpz_clear(n_factors);
    fmpz_clear(f);
}

static void
earlycleanup2(FILE *G, FILE *H, FILE *F, fmpz *sqrts, fmpz *factors, fmpz_t f, fmpz_t n_primes, fmpz_t n_factors)
{
    fclose(H);
    fclose(F);
    fclose(G);
    fmpz_clear(n_primes);
    fmpz_clear(n_factors);
    fmpz_clear(f);
    for (ulong i = 0; i < fmpz_get_ui(n_primes); i++)
    {
        fmpz_clear(sqrts + i);
    }
    for (ulong i = 0; i < fmpz_get_ui(n_factors); i++)
    {
        fmpz_clear(factors + i);
    }
}