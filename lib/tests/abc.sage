from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_frobenius import *
from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
from argparse import ArgumentParser
from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
from sage.schemes.elliptic_curves.isogeny_small_degree import isogenies_2

parser=ArgumentParser()
parser.add_argument('-n', '--nbprimes', default='3')
parser.add_argument('-v', '--verbose', action='store_true')

args=parser.parse_args()

"""
F=open(prime_name, "r")
p=Integer(F.readline())
"""

L.<i>=NumberField(x^2+1)

"""
Goal of this program, Compute h(D) the class number of Z+f*Z[i], for that 
compute some primitive orientation for E0. Find p such that f-torsion is 
defined over F_p.
"""

l=Integer(1050).next_prime()

L1 = L.primes_above(l)[0]
L2 = L.primes_above(5)[0]
L3=(L1**2)*L2

alpha=(L3).gens_reduced()[0]
(f_re, f_im) = (Integer(alpha.real()), Integer(alpha.imag()))
(f_re_abs, f_im_abs)=(abs(f_re), abs(f_im))

"""
For some big l prime in Z and ideals, the primes above l in Z[i].
Check if the generator ideals[0] (mod +-1, +-i). has imaginary of real
part prime. Then Z[generator] is of conductor the prime. The norm will be l."""
while not (\
	(f_re_abs.is_pseudoprime() and l&3==1)\
		or\
	(f_im_abs.is_pseudoprime() and l&3==1) ):

	l=l.next_prime()
	L1 = L.primes_above(l)[0]
	L3=(L1**2)*L2


	alpha=L3.gens_reduced()[0]
	(f_re, f_im) = (Integer(alpha.real()), Integer(alpha.imag()))
	(f_re_abs, f_im_abs)=(abs(f_re), abs(f_im))


if f_re_abs.is_pseudoprime():
	if f_re>0:
		f=f_re
		a=-f_im
	else:
		f=-f_re
		a=f_im
else:
	if f_im>0:
		f=f_im
		a=f_re
	else:
		f=-f_im
		a=-f_re

norm = a**2+f**2
trace = 2*a

alpha=a+i*f
if args.verbose:
	print(f"alpha={alpha} with factored norm: {factor(norm)}\n")

h_f = f-kronecker(-1, f)
if args.verbose:
	print(f"Class number for Z+f*Z[i] with f={f}: {h_f}\n")

"""Find p such that the f-torsion and l-torsion are defined over F_p"""
c=1
p=c*f*l*5-1
pmone=1
"""So that f and l-torsion are defined over F_p"""
while True:
	if p.is_pseudoprime() and p&3==3:
		pmone=1
		break
	p=p+2
	if p.is_pseudoprime() and p&3==3:
		pmone=-1
		break
	c+=1
	p=c*f*l*5-1
	
if args.verbose:
	print(f"Found prime for base field of the curve: p=c*f*l+-1={c}*{f}*{5}*{l}+-1={p}\n")
	print(f"Factored p-({-pmone}): {factor(p+pmone)}")


"""
Next:
	-find a primitively oriented curve by Z+f*Z[i]
"""

K.<t>=GF((p,2), modulus=[1,0,1])

K_l=GF(l)
R.<X>=K_l[]

E0=EllipticCurve(K, [1,0])
E0.set_order((p+1)**2)


P_5,Q_5=E0.torsion_basis(5)
P_f,Q_f=E0.torsion_basis(f)
P_l,Q_l=E0.torsion_basis(l)

iota=WeierstrassIsomorphism(E0,[-t,0,0,0], E0)


"""
From the (non primitive) orientation by alpha on E0, push it through phi_f for some 
phi_f of degree f to get a primitive orientation by alpha on some E.
"""

"""
if w_0=a+iota*f then w_0 has minimal polynomial, (X-(a+sqrt_l*f))(X-(a-sqrt_l*f)).
So that to compute w_0 we can simply find an eigenvector of iota mod l. 
"""

"""
WARNING! GOT STUCK A LONG TIME BECAUSE THE l or l_bar/ 5 or 5_bar 
combination is important
"""
K1.<i>=NumberField(x^2+1)
sqrt_l=K1.primes_above(l)[0].gens_two()[1].real()
sqrt_5=K1.primes_above(5)[0].gens_two()[1].real()
O=P_l-P_l
(P_l_bar, Q_l_bar)=(sqrt_l*P_l+iota(P_l), sqrt_l*Q_l+iota(Q_l))
(P_l, Q_l)=(sqrt_l*P_l-iota(P_l), sqrt_l*Q_l-iota(Q_l))
(P_5, Q_5)=(sqrt_5*P_5-iota(P_5), sqrt_5*Q_5-iota(Q_5))

"""Now P_l, Q_l\in ker(iota-sqrt_l), they even generate it if nonzero (one of them must be)."""

if P_l==O:
	P_l=Q_l

if P_l_bar==O:
	P_l_bar=Q_l_bar

"""Same for P_5, Q_5"""
if P_5==O:
	P_5=Q_5




"""
Compute any isogeny of degree f (if it is not the one computed with the orientation 
we are good (shouldn't happen). To get that the orientation on E_f is the pushforward 
of the one on E0, f should be prime with n(alpha) which is 5*l^2 (and f!=5,l).
"""

phi_f=E0.isogeny(P_f)
E_f=phi_f.codomain()
j=E_f.j_invariant()

if args.verbose:
	print(f"On elliptic curve with j invariant {j}\n")

phi_l=E_f.isogeny(phi_f(P_l))
phi_l_bar=E_f.isogeny(phi_f(P_l_bar))


E_5=phi_l.codomain()
E_l_bar=phi_l_bar.codomain()
phi_5=E_5.isogeny(phi_l(phi_f(P_5)))

E_5_bis=phi_5.codomain()
iso=E_5_bis.isomorphism_to(E_l_bar)

phi_l_bar=phi_l_bar.dual()
w=phi_l_bar*iso*phi_5*phi_l

print(w)




