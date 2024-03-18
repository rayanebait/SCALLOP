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

K.<t>=GF((p,2), modulus=[1,0,1])
R.<x>=K[]

E0=EllipticCurve(K, [1,0])
E0.set_order((p+1)**2)


#print(psi.rational_maps())
#print(phi.rational_maps())


iota=WeierstrassIsomorphism(E0,[-t,0,0,0], E0)

#E0=E0.change_ring(F)
#frob=EllipticCurveHom_frobenius(E0)
#E0=E0.change_ring(K)

#print(frob+iota)

#print(Q, "maps to: ",(frob*iota)(Q))

n=10
m=11

psi_l=E0.division_polynomial(13,x)
psi_iota_l=psi_l(-x)

print(psi_l-psi_iota_l, 5*gcd(psi_l, psi_iota_l))

