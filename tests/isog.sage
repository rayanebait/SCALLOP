from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_frobenius import *
from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
#from sage.schemes.elliptic_curves.ell_finite_field import *
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-n', '--nbprimes', default='14')

args=parser.parse_args()

n=args.nbprimes
prime_name="../txt/prime_"+n+"_primes.md"

n=Integer(n)

F=open(prime_name, "r")
p=Integer(F.readline())

K.<i>=GF((p,2), modulus=[1,0,1])

E0=EllipticCurve(K, [1,0])
E0.set_order((p+1)**2)

print(E0.is_supersingular())

Q=E0.lift_x(4)
P=E0.lift_x(0)

isom=E0.isomorphism_to(E0)
print(isom)

psi=EllipticCurveIsogeny(E0,Q)
phi=EllipticCurveIsogeny(E0,P)

#print(psi.rational_maps())
#print(phi.rational_maps())


iota=WeierstrassIsomorphism(E0,[-i,0,0,0], E0)

#E0=E0.change_ring(F)
frob=EllipticCurveHom_frobenius(E0)
#E0=E0.change_ring(K)

#print(frob+iota)

print(Q, "maps to: ",(frob*iota)(Q))

n=10
m=11

n_=E0.scalar_multiplication(n)
m_=E0.multiplication_by_m(m, x_only=True)

print(m_[0])

