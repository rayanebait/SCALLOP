from sage.schemes.elliptic_curves.weierstrass_morphism import *
from argparse import ArgumentParser
from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius

class EvalIsogSum:
	def __init__(self,precurve, postcurve, basis=[], coeffs=[1,1,1,1]):
		self.precurve=precurve
		self.postcurve=postcurve
		self.basis=basis
		self.coeffs=coeffs

	def fill_basis(self, basis):
		self.basis=basis
	def fill_coeffs(self, coeffs):
		self.coeffs=coeffs
	def evaluate(self, P):
		return sum([self.coeffs[i]*self.basis[i](P) for i in range(len(self.basis))])

proof.arithmetic(False)

parser = ArgumentParser()

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-n', '--nbprimes', default='14')

args = parser.parse_args()
n=args.nbprimes
if n!='3' and n!='14' and n!='19' and n!='20' and n!='26':
	raise SystemExit(f"Currently supports 14, 19, 20 and 26 primes parameters, not {n}\n")

nb_primes=Integer(n)

if args.verbose:
	print(f"Generating oriented curve for {nb_primes} primes parameter\n")

endo_file="../../txt/endo_"+str(nb_primes)+"_primes.md"
prime_file="../../txt/prime_"+str(nb_primes)+"_primes.md"
cond_file="../../txt/conductor_"+str(nb_primes)+"_primes.md"
alpha_file="../../txt/alpha_"+str(nb_primes)+"_primes.md"

F=open(endo_file, "r")
G=open(prime_file, "r")
H=open(cond_file, "r")
O=open(alpha_file, "r")

p=Integer(G.readline())
K.<i>=GF((p,2), modulus=[1,0,1])
R.<x>=K[]

conductor=Integer(O.readline())
alpha2=-conductor
alpha1=Integer(O.readline())

norm=alpha1**2+alpha2**2
trace=2*alpha1


G.close()
H.close()
O.close()

li=5
primes=[]
i=0
while i<nb_primes:
	if li&3==1:
		primes.append(li)
		i+=1
	li=li.next_prime()


E0=EllipticCurve(K, [1,0])
pols=[]
psi_alpha1_l=1
psi_alpha2_l=1
alpha1_l=1
alpha2_l=1

i=0
for l in primes:
	if i>5:
		break
	i+=1

	alpha1_l=Mod(alpha1,l)
	alpha2_l=Mod(alpha2,l)
	
	#gcd(psi_alpha1_l(x),psi_alpha2_l(-x), x^(p**2)-x) ?
	psi_alpha1_l=E0.division_polynomial(alpha1_l)
	psi_alpha2_l=E0.division_polynomial(alpha2_l)


	xbar=R.quotient_by_principal_ideal((psi_alpha1_l)*R).0
	pols.append((psi_alpha1_l, psi_alpha2_l(-x), (xbar**(p)-x) ))

print( gcd(gcd(pols[1][0], pols[1][1]), pols[1][2]) )
F.close()


#iota=WeierstrassIsomorphism(E0, [-i,0,0,0], E0)

#P=E0.lift_x(4)
#Q=iota(P)
#E0.division_polynomial(10,x)



