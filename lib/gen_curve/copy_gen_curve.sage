from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite hom_comp
from argparse import ArgumentParser


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


	
#Shouldn't be heavy
#def compute_kernel_sum(a1_l, a2_l, l, E0):

	#if a1_l==0:
		#return E0.division_polynomial(a2_l)(-x)
	#if a2_l==0:
		#return E0.division_polynomial(a1_l)

	#ADD a1_l==a2_l CASE

	#psi_a1_l=[R(E0.division_polynomial(a1_l+i, x)) if a1_l+i!=0 else R(0) for i in IntegerRange(-2,3)]
	#psi_a2_l=[R(E0.division_polynomial(a2_l+i, x)(-x)) if a2_l+i!=0 else R(0) for i in IntegerRange(-2,3)]

	#Mult by a1_l
	#x_1=x-(psi_a1_l[1]*psi_a1_l[3])/psi_a1_l[2]**2
	#y_1=(psi_a1_l[4]*psi_a1_l[1]**2-psi_a1_l[0]*psi_a1_l[3]**2)/(4*psi_a1_l[2]**3)

	#Mult by a2_l*iota
	#x_2=-x-(psi_a2_l[1]*psi_a2_l[3])/psi_a2_l[2]**2
	#y_2=(psi_a2_l[4]*psi_a2_l[1]**2-psi_a2_l[0]*psi_a2_l[3]**2)/(4*psi_a2_l[2]**3)

	#lambd=((y_2-y_1)/(x_2-x_1))**2
	
	#x_sum=-x_1-x_2+lambd
	#y_sum=-y_1-lambd*(x_sum-x_1)

	#left_x_summand=(psi_a2_l[2]**2)*(x*(psi_a1_l[2]**2)-psi_a1_l[1]*psi_a1_l[3])
	#right_x_summand=(psi_a1_l[2]**2)*(-x*(psi_a2_l[2]**2)-psi_a2_l[1]*psi_a2_l[3])

	#left_y_summand=(psi_a1_l[4]*psi_a1_l[1]**2-psi_a1_l[0]*psi_a1_l[3]**2)*(4*psi_a2_l[2]**3)
	#right_y_summand=(psi_a2_l[4]*psi_a2_l[1]**2-psi_a2_l[0]*psi_a2_l[3]**2)*(4*psi_a1_l[2]**3)

	#return left_x_summand-right_x_summand

O=E0(0,0)-E0(0,0)
print(factor(norm))



def compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm):
	i=0
	gen_of_l_part=[]
	isogeny_chain=[]
	curve_chain=[]
	Ei=E0
	for l in primes:
		print(f"Prime:{l}\n")
		if i>3:
			break
		i+=1
	
		#psi_l=E0.division_polynomial(l)
		P, Q=E0.torsion_basis(l)
		sqrt_l=square_root_mod_prime(Mod(-1, l),l)
		if sqrt_l==0:
			a=E0.scalar_multiplication(alpha1-trace)
			b=E0.scalar_multiplication(alpha2)
		else :
			a=E0.scalar_multiplication(alpha1-floor(norm/sqrt_l))
			b=E0.scalar_multiplication(alpha2)

		Pi=a(P)+iota*b(P)
		if (a(P_)+iota*b(P_))==O:
			gens_of_l_part.append(Pi)
		else:
			Qi=a(Q)+iota*b(Q)
			if (a(Qi)+iota*b(Qi))==O:
				gens_of_l_part.append(Qi)

	for i in range(len(primes)):
		P=gen_of_l_part[i]
		for phi in isogeny_chain:
			P=phi(P)

		phi_i=EllipticCurveIsogeny(curve_chain[i],P)
		curve_chain.append(phi_i.codomain())
		isogeny_chain.append(phi_i)

	phi_L1L2=hom_comp.from_factors(isogeny_chain, E0)

	

	#Rbar.<xbar>=R.quotientt(psi_l)
	#frob_l=Rbar.lift(xbar**(p*p))-x
	#print(frob_l)

	#ker=compute_kernel_sum(alpha1_l, alpha2_l, l, E0)
	#print("l-torsion part in the kernel: ", gcd(ker, psi_l))
	
F.close()
