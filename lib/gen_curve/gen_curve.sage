from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
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
#(x,y)->(-x, iy)
iota=WeierstrassIsomorphism(E0, [-J,0,0,0], E0)


	

print(factor(norm))
print(p)
print(f"p%{5}={p%5}, p%{13}={p%13}\n")



#Compute the kernel of w_0 by writing it as phi_L1L2*phi_L1^-1_dual and 
#compute phi_L1 as the chain of phi_li_Ei. Compute an 
#the kernel of phi_li by computing the kernel of w_0-sqrt_l 
#as done in SCALLOP below proposition 8 page 11
def compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm):
	O=E0(0,0)-E0(0,0)

	gen_of_l_part=[]
	gen_of_l_part_bar=[]
	isogeny_chain=[]
	isogeny_chain_bar=[]
	curve_chain=[]
	curve_chain_bar=[]

	#Compute a generator of the kernel of each E[I_l]
	for l in primes:
		print(f"Prime:{l}\n")
		print(E0.division_field(l))
	
		r=0
		P, Q=E0.torsion_basis(l)
		sqrt_l=square_root_mod_prime(Mod(-1, l),p=l)

		print(sqrt_l)
		a_l=E0.scalar_multiplication(alpha1-sqrt_l)
		b=E0.scalar_multiplication(alpha2)
		bbar=E0.scalar_multiplication(-alpha2)
		if sqrt_l==0:
			a=E0.scalar_multiplication(alpha1-trace)
		else:
			#plutot multiplier par l'inverse de sqrt_l mod l ?
			#w_0=sqrt_l sur E[I] i.e tr(w_0)=sqrt_l+n(alpha)*sqrt_l^-1 mod l?

			#inv=ZZ(GF(l)(sqrt_l)**(-1))
			correc=Mod(trace-sqrt_l, l)
			inv=floor(norm/sqrt_l)
			a=E0.scalar_multiplication(alpha1-correc)

		P_=a(P)+(iota*b)(P)
		if (a_l(P_)+(iota*b)(P_))!=O:
			r+=1
			Q_=a(Q)+(iota*b)(Q)
			if (a_l(Q_)+(iota*b)(Q_))!=O:
				raise RuntimeError(f"Curve Generation failed: Couldn't find generator for E[L1L2], round {r}, prime {l}\n")
			gen_of_l_part.append(Q_)
		else:
			gen_of_l_part.append(P_)

		P_bar=a(P)+(iota*bbar)(P)
		if (a_l(P_bar)+(iota*bbar)(P_bar))!=O:
			Q_bar=a(Q)+(iota*bbar)(Q)
			if (a_l(Q_bar)+(iota*bbar)(Q_bar))!=O:
				raise RuntimeError("Curve Generation failed: Couldn't find generator for E[L1]\n")
			gen_of_l_part_bar.append(Q_bar)
		else:
			gen_of_l_part_bar.append(P_bar)


	phi_5=EllipticCurveIsogeny(E0,gen_of_l_part[0])
	print(phi_5)
	curve_chain.append(phi_5.codomain())
	isogeny_chain.append(phi_5)
	curve_chain_bar.append(E0)

	#compute the isogeny factors of phi_L1L2 and phi_L1^-1
	for i in range(1,len(primes)):
		P=gen_of_l_part[i]
		Pbar=gen_of_l_part_bar[i]
		#push a generator of E0[I_l] through every isogeny to get E_i[I_l]
		for (phi,phibar) in zip(isogeny_chain, isogeny_chain_bar):
			P=phi(P)
			Pbar=phibar(Pbar)

		phi_i=EllipticCurveIsogeny(curve_chain[i-1],P)
		phi_i_bar=EllipticCurveIsogeny(curve_chain_bar[i-1],Pbar)
		print(phi_i)

		curve_chain.append(phi_i.codomain())
		curve_chain_bar.append(phi_i_bar.codomain())
		isogeny_chain.append(phi_i)
		isogeny_chain_bar.append(phi_i_bar)

	phi_L1L2=hom_comp.from_factors(isogeny_chain, E0)
	phi_L1bar=(hom_comp.from_factors(isogeny_chain_bar, E0)).dual()

	w_0=phi_L1bar*phi_L1L2

	print(factor(phi_L1L2.degree()))

	

compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm)



F.close()
