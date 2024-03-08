#ifndef qfb_H
#define qfb_H

#define qfb_INLINE static __inline__

#include "fmpz.h"


typedef struct{
    fmpz_t a;
    fmpz_t b;
    fmpz_t c;
} qfb;

typedef qfb qfb_t[1];

qfb_INLINE
void qfb_init(qfb_t q)
{
   fmpz_init(q->a);
   fmpz_init(q->b);
   fmpz_init(q->c);
}

qfb_INLINE
void qfb_clear(qfb_t q)
{
   fmpz_clear(q->a);
   fmpz_clear(q->b);
   fmpz_clear(q->c);
}

qfb_INLINE
void qfb_setcoef1(qfb_t q, fmpz_t a)
{
   fmpz_set(q->a, a);
}

qfb_INLINE
void qfb_setcoef2(qfb_t q, fmpz_t b)
{
   fmpz_set(q->b, b);
}

qfb_INLINE
void qfb_setcoef3(qfb_t q, fmpz_t c)
{
   fmpz_set(q->c, c);
}

qfb_INLINE
void qfb_setcoeffs(qfb_t q, fmpz_t a, fmpz_t b, fmpz_t c)
{
   fmpz_set(q->a, a);
   fmpz_set(q->b, b);
   fmpz_set(q->c, c);
}

qfb_INLINE
int qfb_equal(qfb_t f, qfb_t g)
{
   return (fmpz_equal(f->a, g->a)
        && fmpz_equal(f->b, g->b)
        && fmpz_equal(f->c, g->c));
}

qfb_INLINE
void qfb_set(qfb_t f, qfb_t g)
{
   fmpz_set(f->a, g->a);
   fmpz_set(f->b, g->b);
   fmpz_set(f->c, g->c);
}

qfb_INLINE
void qfb_discriminant(fmpz_t D, qfb_t f)
{
   fmpz_t t;
   fmpz_init(t);

   fmpz_mul(t, f->a, f->c);
   fmpz_mul_2exp(t, t, 2);
   fmpz_mul(D, f->b, f->b);
   fmpz_sub(D, D, t);

   fmpz_clear(t);
}

qfb_INLINE
void qfb_print(qfb_t q)
{
    printf("(");
    fmpz_print(q->a);
    printf(", ");
    fmpz_print(q->b);
    printf(", ");
    fmpz_print(q->c);
    printf(")");
}

qfb_INLINE
void qfb_array_clear(qfb **forms, slong num)
{
   slong k;

   for (k = 0; k < num; k++)
   {
      fmpz_clear((*forms)[k].a);
      fmpz_clear((*forms)[k].b);
      fmpz_clear((*forms)[k].c);
   }
   flint_free(*forms);
}


void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D);

qfb_INLINE
int qfb_is_reduced(qfb_t q)
{
   fmpz_t tmp;
   fmpz_init(tmp);
   fmpz_abs(tmp, q->b);

   if(fmpz_cmp(tmp, q->a)<=0 && fmpz_cmp(q->a, q->c)<=0){
      if(fmpz_cmp(tmp, q->a) == 0 || fmpz_cmp(q->a, q->c) == 0){
         if(fmpz_cmp_ui(q->b, 0) >= 0){
            fmpz_clear(tmp);
            return 1;
         }
      }
   }

   fmpz_clear(tmp);
   return 0;
}

/*slong qfb_reduced_forms(qfb ** forms, slong d);

slong qfb_reduced_forms_large(qfb ** forms, slong d);*/

void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t D, fmpz_t L);

void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L);

void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp);

void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp);

void qfb_pow_with_root(qfb_t r, qfb_t f, fmpz_t D, fmpz_t e, fmpz_t L);

qfb_INLINE
void qfb_inverse(qfb_t r, qfb_t f)
{
   qfb_set(r, f);

   if (fmpz_equal(f->a, f->c)
    || fmpz_equal(f->a, f->b))
    return;

   fmpz_neg(r->b, r->b);
}

qfb_INLINE
int qfb_is_principal_form(qfb_t f, fmpz_t D)
{
   if (!fmpz_is_one(f->a))
      return 0;

   if (fmpz_is_odd(D)) /* D = 1 mod 4 */
      return fmpz_is_one(f->b);

   return fmpz_is_zero(f->b); /* D = 0 mod 4 */
}

qfb_INLINE
void qfb_principal_form(qfb_t f, fmpz_t D)
{
   fmpz_set_ui(f->a, 1);

   if (fmpz_is_odd(D)) /* D = 1 mod 4 */
      fmpz_set_ui(f->b, 1);
   else /* D = 0 mod 4 */
      fmpz_set_ui(f->b, 0);

   fmpz_sub(f->c, f->b, D);
   fmpz_fdiv_q_2exp(f->c, f->c, 2);
}

qfb_INLINE
int qfb_is_primitive(qfb_t f)
{
   fmpz_t g;
   int res;

   fmpz_init(g);
   fmpz_gcd3(g, f->a, f->b, f->c);
   res = fmpz_is_pm1(g);
   fmpz_clear(g);

   return res;
}

void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p);

ulong find_power(qfb_t f, fmpz_t n, ulong base);

int qfb_exponent_element(fmpz_t exponent, qfb_t f, fmpz_t n, ulong B1, ulong B2_sqrt);

#endif
