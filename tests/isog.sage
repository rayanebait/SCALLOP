from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_frobenius import *
from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
from argparse import ArgumentParser
from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
from sage.schemes.elliptic_curves.isogeny_small_degree import isogenies_2

parser=ArgumentParser()
parser.add_argument('-n', '--nbprimes', default='3')

args=parser.parse_args()



n=args.nbprimes
prime_name="../txt/prime_"+n+"_primes.md"

n=Integer(n)

F=open(prime_name, "r")
p=Integer(F.readline())

K.<t>=GF((p,2), modulus=[1,0,1])
R.<x>=K[]

E0=EllipticCurve(K, [1,0])
E0.set_order((p+1)**2)



iota=WeierstrassIsomorphism(E0,[-t,0,0,0], E0)
l=13
P,Q=E0.torsion_basis(l)

sqrt_l=square_root_mod_prime(Mod(-1, l), p=l)

mul_sqrt_l=E0.scalar_multiplication(sqrt_l)

R=E0(2627*t + 3170, 8833*t + 4122, 1)

print(f"{R}")
print(f"{mul_sqrt_l(P)}, {iota(P)}")
print(f"{mul_sqrt_l(Q)}, {iota(Q)}")
print(f"{mul_sqrt_l(Q)==iota(Q)}")
print(f"{P.order()}, {Q.order()}")


