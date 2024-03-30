from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
from sage.schemes.elliptic_curves.isogeny_small_degree import isogenies_2
from argparse import ArgumentParser


proof.arithmetic(False)

parser = ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-n', '--nbprimes', default='3')
parser.add_argument('-s', '--samples', default='100')
args = parser.parse_args()


n=args.nbprimes
if n!='3' and n!='14' and n!='19' and n!='20' and n!='26':
	raise SystemExit(f"Currently supports 3, 14, 19, 20 and 26 primes parameters, not {n}\n")

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


"""
	Compute the kernel of w_0 by writing it as phi_L1L2*phi_L1^-1_dual and 
	compute phi_L1 as the chain of phi_li_Ei. Compute an 
	the kernel of phi_li by computing the kernel of w_0-sqrt_l 
	as done in SCALLOP below proposition 8 page 11.
"""
def compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm):
	O=E0(0,0)-E0(0,0)

	gen_of_l_part=[]
	gen_of_l_part_bar=[]
	isogeny_chain=[]
	isogeny_chain_bar=[]
	curve_chain=[]
	curve_chain_bar=[]

	#Compute a generator of the kernel of each E[I_l]
	if args.verbose:
		print(f"Computing orientation w_0=alpha1+i*alpha2 with \n\
			\talpha1={alpha1}\n\
			\talpha2={alpha2}\n\
			And norm: {factor(norm)}")

	for l in primes[3:]:
		if args.verbose:
			print(f"\nComputing generator for E[I] for I a prime above prime :{l}\n")
	
		r=0
		P, Q=E0.torsion_basis(l)
		if args.verbose:
			l_mul=E0.scalar_multiplication(l)
			print(f"Generators for the l-torsion:\n\t{P}, {Q}\n\nWith [l]*{P}={l_mul(P)}\n\nand [l]*{Q}={l_mul(Q)}\n")
		sqrt_l=square_root_mod_prime(Mod(-1,l),p=l)


		a_l=E0.scalar_multiplication(Mod(alpha1, l)-sqrt_l)
		b=E0.scalar_multiplication(Mod(alpha2,l))
		bbar=E0.scalar_multiplication(-Mod(alpha2, l))
		if sqrt_l==0:
			a=E0.scalar_multiplication(Mod(alpha1-trace, l))
		else:
			#inv=ZZ(GF(l)(sqrt_l)**(-1))
			correc=Mod(trace-sqrt_l, l)
			inv=norm/sqrt_l
			a=E0.scalar_multiplication(Mod(alpha1-inv, l))
			print(f"norm : {norm}")
			print(f"Trace mod l: {Mod(trace,l)}\n +-sqrt_l: {sqrt_l}, {Mod(-sqrt_l,l)}\n +-inv: {inv}, {Mod(-inv,l)}\n")
			return

		if args.verbose:
			print(f"Pushing {P} through alpha1+ialpha2")
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

	kernel_iso_components=[]
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

	#should only compute the kernel 
	phi_L1L2=hom_comp.from_factors(isogeny_chain, E0)
	phi_L1bar=(hom_comp.from_factors(isogeny_chain_bar, E0)).dual()

	w_0=phi_L1bar*phi_L1L2

	print(factor(phi_L1L2.degree()))

	
endo=[]
for i in range(4):
	endo.append(Integer(F.readline()))

l_0=Integer(F.readline())
h=Integer(F.readline())

#only done for l_0=2
#very slow
"""
E0->isogs[0].codomain()->isogs[0].codomain()->....
Pour le parcours en entier:
y'a une bijection -> chemin d'un arbre 3-régulier de profondeur h <-> {0,1,2}^h
(y'a de la redondance à cause des multiplications mais pas possible de savoir 
si isogs[0]: E1->E2 = isogenies_2(isogs[0].codomain())[i].dual() sans faire le calcul)

Pour la fonction next:
On a un entier en base 3, a_0...a_h si on est à l'étage i on avance jusqu'au bout de l'entier,
puis on recule et tant que le digit est =2 on incrémente l'entier et on recule, si on a pas deux,
on suit le chemin jusqu'au bout et on recommence. En mémoire faut garder une map E_i->2-isogénies
sortantes. Et le chemin entier et l'ensemble des isogénies de degré 2^h calculé .
"""

"""
For the computation of l_0^h isogenies:
	Visualise paths of length h in the l_0+1 regular graph as 
	integers in base l_0+1. Integers in base l_0+1 are in bijection
	with homomorphism from l_0+1 regular trees of depth h so that 
	we may get the next path as:
		-go to depth h
		-while a_h!=l_0: a_h+=1 (visit all l_0 leaves)
		-while a_i==l_0: a_i=0, i-=1 (go to the next subtree)
	which in base l_0+1 means just incrementing the integer
"""



def brute_force_orientation_small(E0, primes, alpha1,\
				alpha2, trace, norm,\
				l_0, h, l0_h_isogs=[],\
				path=[], path_base_l0_plus_1=[],\
				step=0, E_to_l0_isogs={},\
				intermediate_stop=100):

	def go_to_depth(h, path,\
		path_base_l0_plus_1,\
		step, E_to_l0_isogs,\
		E=E0, init=False):

		while step < h-1:
			"""
			Check if the codomain has already been encountered
			"""
			l0_isogs=E_to_l0_isogs.get(E)
			if l0_isogs == None:
				l0_isogs=isogenies_2(E)
				E_to_l0_isogs[E]=l0_isogs
			if init:
				path.append(l0_isogs[path_base_l0_plus_1[step]])
			else: 
				path[step]=l0_isogs[path_base_l0_plus_1[step]]

			E=path[step].codomain()
			step+=1

	def add_isogs(base, path,\
		E, l0_h_isogs,\
		E_to_l0_isogs,\
		nb_isogs):

		l0_isogs=E_to_l0_isogs.get(E)
		if l0_isogs==None:
			l0_isogs=isogenies_2(E)
			E_to_l0_isogs[E]=l0_isogs
		

		for i in range(base):
			l0_h_isogs.append(hom_comp.from_factors(path+[l0_isogs[i]]))
		nb_isogs+=base
		
		
	def go_back(base, h, step,\
		path_base_l0_plus_1,
		intermediate_stop=0):

		while path_base_l0_plus_1[step]==base-1 and step>=intermediate_stop:
			path_base_l0_plus_1[step]=0
			step-=1
		if step==intermediate_stop:
			return True
		else:
			return False

	E=E0
	go_to_depth(h, path,\
		path_base_l0_plus_1,\
		step, E_to_l0_isogs, init=True)

	done=False

	nb_isogs=0
	while not done:
		go_to_depth(h, path,\
			path_base_l0_plus_1,\
			step, E_to_l0_isogs,\
			E=E)

		add_isogs(l_0+1,path,\
			E, l0_h_isogs,\
			E_to_l0_isogs,\
			nb_isogs)
		
		done=go_back(h, step,\
			path_base_l0_plus_1,\
			intermediate_stop)

	return isogs


	
	

def brute_force_orientation(E0, primes, alpha1, alpha2, trace, norm):
	G=cartesian_product([[0,1,2]]*h)
	for seq in G:
		isog_factors=[]
		E_j=E0
	#should find a way to not recompute every l_0 isogs each time
		i=0
		for j in seq:
			isogs=isogenies_2(E_j)
			isog_factors.append(isogs[j])
			E_j=isogs[j].codomain()
		#only compute kernel
		Phi=hom_comp.from_factors(isog_factors, E0)	
		print(Phi)
		print(Phi.degree())

	

#compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm)

nb_samples=Integer(args.samples)
isogs=brute_force_orientation_small(E0, primes, alpha1,\
				alpha2, trace, norm,\
				l_0, h, path_base_l0_plus_1=[0]*h,\
				intermediate_stop=nb_samples)

print(isogs)

F.close()
