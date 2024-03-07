#include <qfb.h>
#include <fmpz.h>
#include <fmpz_mat.h>


void
qfb_init_set(qfb_t q, fmpz_t a, fmpz_t b, fmpz_t c){
    fmpz_init_set(q->a, a);
    fmpz_init_set(q->b, b);
    fmpz_init_set(q->c, c);
}


/*Only for I prime in Z[i] with I=[1+r^2, f(i-r)] and r^2=-1 mod N(I)*/
void
ideal_to_qfb(qfb_t q, fmpz_t r, fmpz_t f){
    fmpz_t a,b,c;
    fmpz_init_ui(a, 1);
    fmpz_init(b, r);
    fmpz_init(c, f);

    fmpz_addmul(a, b, b);
    fmpz_mul(b, b, c);
    fmpz_mul_ui(b, b, 2);

    fmpz_mul(c, f,f);

    qfb_init_set(q, a,b,c);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}


/*Currently uses rho-Pollard*/
qfb_discrete_log_prime_order(fmpz_t pdlog, qfb_t h, qfb_t g, fmpz_t p)
{

}

/*Computes Pohlig hellman reduction, factors should be ordered*/
void
qfb_discrete_log(fmpz_t dlog, qfb_t h, qfb_t g, fmpz_t f, fmpz *factors, ulong nb_factors)
{
    fmpz_t g_i;
    fmpz_init_set(g_i, factors);
    for(ulong i = 1; i<nb_factors; i++){
    }

}


