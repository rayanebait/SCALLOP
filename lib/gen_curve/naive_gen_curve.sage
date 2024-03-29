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
		return sum([self.coeffs[i]*self.basis[j](P) for j in range(len(self.basis))])


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
K.<J>=GF((p,2), modulus=[1,0,1])
R.<x>=K[]
FracR=FractionField(R)

conductor=Integer(O.readline())
alpha2=conductor
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
E0.set_order((p+1)**2)
iota=WeierstrassIsomorphism(E0, [-J,0,0,0], E0)

alpha1_l=1
alpha2_l=1


	
#Shouldn't be heavy
def compute_kernel_sum(a1_l, a2_l, l, E0):

	if a1_l==0:
		return E0.division_polynomial(a2_l)(-x)
	if a2_l==0:
		return E0.division_polynomial(a1_l)

	#ADD a1_l==a2_l CASE

	psi_a1_l=[R(E0.division_polynomial(a1_l+i, x)) if a1_l+i!=0 else R(0) for i in IntegerRange(-2,3)]
	psi_a2_l=[R(E0.division_polynomial(a2_l+i, x)(-x)) if a2_l+i!=0 else R(0) for i in IntegerRange(-2,3)]

	#Mult by a1_l
	#x_1=x-(psi_a1_l[1]*psi_a1_l[3])/psi_a1_l[2]**2
	#y_1=(psi_a1_l[4]*psi_a1_l[1]**2-psi_a1_l[0]*psi_a1_l[3]**2)/(4*psi_a1_l[2]**3)

	#Mult by a2_l*iota
	#x_2=-x-(psi_a2_l[1]*psi_a2_l[3])/psi_a2_l[2]**2
	#y_2=(psi_a2_l[4]*psi_a2_l[1]**2-psi_a2_l[0]*psi_a2_l[3]**2)/(4*psi_a2_l[2]**3)

	#lambd=((y_2-y_1)/(x_2-x_1))**2
	
	#x_sum=-x_1-x_2+lambd
	#y_sum=-y_1-lambd*(x_sum-x_1)

	left_x_summand=(psi_a2_l[2]**2)*(x*(psi_a1_l[2]**2)-psi_a1_l[1]*psi_a1_l[3])
	right_x_summand=(psi_a1_l[2]**2)*(-x*(psi_a2_l[2]**2)-psi_a2_l[1]*psi_a2_l[3])

	#left_y_summand=(psi_a1_l[4]*psi_a1_l[1]**2-psi_a1_l[0]*psi_a1_l[3]**2)*(4*psi_a2_l[2]**3)
	#right_y_summand=(psi_a2_l[4]*psi_a2_l[1]**2-psi_a2_l[0]*psi_a2_l[3]**2)*(4*psi_a1_l[2]**3)

	return left_x_summand-right_x_summand


print(factor(norm))

i=0
for l in primes:
	print(f"Prime:{l}\n")
	if i>3:
		break
	i+=1

	alpha1_l=Mod(alpha1,l)
	alpha2_l=Mod(alpha2,l)
	psi_l=E0.division_polynomial(l)

	#Rbar.<xbar>=R.quotientt(psi_l)
	#frob_l=Rbar.lift(xbar**(p*p))-x
	#print(frob_l)

	ker=compute_kernel_sum(alpha1_l, alpha2_l, l, E0)
	#print(psi_l)
	print("l-torsion part in the kernel: ", gcd(ker, psi_l))

	#phi=EllipticCurveIsogeny(E0, ker*(ker.coefficients()[-1]**(-1)))

	#roots=ker.roots()
	#a1_l_mul=E0.scalar_multiplication(alpha1_l)
	#a2_l_mul=E0.scalar_multiplication(alpha2_l)

	#for root in roots:
		#try:
			#P=E0.lift_x(root[0])
			#print(phi(P))
			#print(a1_l_mul(P)[0]==a2_l_mul(iota(P)))

		#except ValueError:
			#continue
	

F.close()




#P=E0.lift_x(4)
#Q=iota(P)
#E0.division_polynomial(10,x)



