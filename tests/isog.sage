from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_frobenius import *
from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
#from sage.schemes.elliptic_curves.ell_finite_field import *

F=GF(95677801353807797939337736858276469382798403327077326287693919441924240579)
K.<t>=GF((95677801353807797939337736858276469382798403327077326287693919441924240579,2), modulus=[1,0,1])

E0=EllipticCurve(K, [1,0])
E0.set_order((95677801353807797939337736858276469382798403327077326287693919441924240579+1)**2)

print(E0.is_supersingular())

Q=E0.lift_x(4)
P=E0.lift_x(0)

isom=E0.isomorphism_to(E0)
print(isom)

psi=EllipticCurveIsogeny(E0,Q)
phi=EllipticCurveIsogeny(E0,P)

#print(psi.rational_maps())
#print(phi.rational_maps())


iota=WeierstrassIsomorphism(E0,[-t,0,0,0], E0)

#E0=E0.change_ring(F)
frob=EllipticCurveHom_frobenius(E0)
#E0=E0.change_ring(K)

#print(frob+iota)

print(Q, "maps to: ",(frob*iota)(Q))

n=10
m=11

n_=E0.scalar_multiplication(n)
m_=E0.scalar_multiplication(m)

print(n_.rational_maps()[0].numerator().category())
