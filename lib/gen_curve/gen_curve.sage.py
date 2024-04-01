

# This file was *autogenerated* from the file gen_curve.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_5 = Integer(5); _sage_const_3 = Integer(3); _sage_const_1728 = Integer(1728); _sage_const_4 = Integer(4)
from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
from sage.schemes.elliptic_curves.isogeny_small_degree import isogenies_2
from argparse import ArgumentParser
import time 


proof.arithmetic(False)

parser = ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-n', '--nbprimes', default='3')
parser.add_argument('-s', '--samples', default='0')
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
K = GF((p,_sage_const_2 ), modulus=[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ], names=('J',)); (J,) = K._first_ngens(1)
R = K['x']; (x,) = R._first_ngens(1)
FracR=FractionField(R)

conductor=Integer(O.readline())

alpha2=conductor
alpha1=Integer(O.readline())

norm=alpha1**_sage_const_2 +alpha2**_sage_const_2 
trace=_sage_const_2 *alpha1


G.close()
H.close()
O.close()

li=_sage_const_5 
primes=[]
i=_sage_const_0 
while i<nb_primes:
	if li&_sage_const_3 ==_sage_const_1 :
		primes.append(li)
		i+=_sage_const_1 
	li=li.next_prime()


E0=EllipticCurve(K, [_sage_const_1 ,_sage_const_0 ])
#(x,y)->(-x, iy)
iota=WeierstrassIsomorphism(E0, [-J,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ], E0)




print(factor(norm))


"""
	Compute the kernel of w_0 by writing it as phi_L1L2*phi_L1^-1_dual and 
	compute phi_L1 as the chain of phi_li_Ei. Compute an 
	the kernel of phi_li by computing the kernel of w_0-sqrt_l 
	as done in SCALLOP below proposition 8 page 11.
"""
def compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm):
	O=E0(_sage_const_0 ,_sage_const_0 )-E0(_sage_const_0 ,_sage_const_0 )

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

	for l in primes:
		if args.verbose:
			print(f"\nComputing generator for E[I] for I a prime above prime :{l}\n")
	
		r=_sage_const_0 
		P, Q=E0.torsion_basis(l)
		if args.verbose:
			print(f"Generators for the l-torsion:\n\t{P}, {Q}\n\n")

		a_l=E0.scalar_multiplication(Mod(alpha1-trace, l))
		a=E0.scalar_multiplication(Mod(alpha1, l))
		b=E0.scalar_multiplication(Mod(alpha2,l))
		bbar=E0.scalar_multiplication(-Mod(alpha2, l))

		if args.verbose:
			print(f"Pushing {P} through w\n")

		P_=a_l(P)+(iota*b)(P)
		if P_==O or (a(P_)+(iota*b)(P_))!=O:
			Q_=a_l(Q)+(iota*b)(Q)
			if Q_==O or (a(Q_)+(iota*b)(Q_))!=O:
				raise RuntimeError(f"Curve Generation failed: Couldn't find generator for E[L1L2], round {r}, prime {l}\n")
			gen_of_l_part.append(Q_)
		else:
			gen_of_l_part.append(P_)

		"""Case l=5 can be jumped here"""
		P_bar=a_l(P)+(iota*bbar)(P)
		if P_bar==O or (a(P_bar)+(iota*bbar)(P_bar))!=O:
			Q_bar=a_l(Q)+(iota*bbar)(Q)
			if Q_bar==O or (a(Q_bar)+(iota*bbar)(Q_bar))!=O:
				raise RuntimeError("Curve Generation failed: Couldn't find generator for E[L1^-1]\n")
			gen_of_l_part_bar.append(Q_bar)
		else:
			gen_of_l_part_bar.append(P_bar)


	E_i, E_i_bar=(E0,E0)
	j_path=[_sage_const_1728 ]
	j_bar_path=[_sage_const_1728 ]


	#compute the isogeny factors of phi_L1L2 and phi_L1^-1
	for (P,Pbar) in zip(gen_of_l_part[_sage_const_1 :], gen_of_l_part_bar[_sage_const_1 :]):
		#push a generator of E0[I_l] through every isogeny to get E_i[I_l]
		for (phi,phibar) in zip(isogeny_chain, isogeny_chain_bar):
			P=phi(P)
			Pbar=phibar(Pbar)

		phi_i=EllipticCurveIsogeny(E_i,P)
		phi_i_bar=EllipticCurveIsogeny(E_i_bar,Pbar)

		E_i=phi_i.codomain()
		E_i_bar=phi_i_bar.codomain()

		j_path.append(E_i.j_invariant())
		j_bar_path.append(E_i_bar.j_invariant())

		isogeny_chain.append(phi_i)
		isogeny_chain_bar.append(phi_i_bar)

	P=gen_of_l_part[_sage_const_0 ]
	for phi in isogeny_chain:
		P=phi(P)

	phi_5=EllipticCurveIsogeny(E_i, P)
	isogeny_chain.append(phi_5)

	j_path.append(phi_5.codomain().j_invariant())

	#should only compute the kernel 
	phi_L1L2=hom_comp.from_factors(isogeny_chain, E0)
	phi_L1bar=(hom_comp.from_factors(isogeny_chain_bar, E0)).dual()
	#phi_L1bar=(hom_comp.from_factors(isogeny_chain[:-1], E0)).dual()

	if args.verbose:
		path_str='--->'.join(str(j) for j in j_path[:len(j_path)-_sage_const_1 ])
		path_str+="--->"+str(j_path[-_sage_const_1 ])

		print(f"Took path:\n\t")
		print(f"{path_str} for phi_L1L2\n\n")

		path_str='--->'.join(str(j) for j in j_bar_path[:len(j_bar_path)-_sage_const_1 ])
		path_str+="--->"+str(j_bar_path[-_sage_const_1 ])

		print(f"Took path:\n\t")
		print(f"{path_str} for phi_L1^-1\n\n")
	try:
		w_0=phi_L1bar*phi_L1L2
	except TypeError:
		E=phi_L1L2.codomain()
		E_=phi_L1bar.domain()
		iso=E.isomorphism_to(E_)
		print(iso)
		w_0=phi_L1bar*iso*phi_L1L2

		return w_0

	
endo=[]
#norm2=0
for i in range(_sage_const_4 ):
	endo.append(Integer(F.readline())//_sage_const_2 )

#print(f"Factored norm of endomorphism {endo}: {factor(norm2)}\n")
l0=Integer(F.readline())
h=Integer(F.readline())

#only done for l_0=2
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
		-while a_h!=l_0: a_h+=1 (visit all l_0+1 leaves)
		-while a_i==l_0: a_i=0, i-=1 (go back to an unvisited subtree)
	which in base l_0+1 means just incrementing the integer
"""



def brute_force_orientation_small(E0, primes, alpha1,				alpha2, trace, norm,				l0, h,				intermediate_stop=_sage_const_0 ):

	if intermediate_stop>h:
		raise ValueError(f"-i should be < to h\n")
		return

	def format_j_path(j_path):
		return "--->".join(str(j) for j in j_path[:len(j_path)-_sage_const_1 ])			+str(j_path[-_sage_const_1 ])

	def go_to_depth(h, path_base_l0_plus_1,			init=False, testing=False):

		nonlocal path, step,			E_to_l0_isogs, E
		if testing:
			nonlocal j_path

		while step < h-_sage_const_1 :
			if args.verbose:
				print(f"At node {E.j_invariant()} at step {step}\n")
			"""
			Check if the codomain has already been encountered
			"""
			l0_isogs=E_to_l0_isogs.get(E)
			if l0_isogs == None:
				l0_isogs=isogenies_2(E)
				E_to_l0_isogs[E]=l0_isogs
			if init:
				path.append(l0_isogs[path_base_l0_plus_1[step]])
				if testing:
					j_path.append(E.j_invariant())
			else: 
				path[step]=l0_isogs[path_base_l0_plus_1[step]]
				if testing:
					j_path[step]=E.j_invariant()
			
			E=path[step].codomain()
			step+=_sage_const_1 

	def add_isogs(base, path, testing=False):
		nonlocal E, E_to_l0_isogs, nb_isogs, j_path
		if testing:
			nonlocal j_path

		l0_isogs=E_to_l0_isogs.get(E)
		if l0_isogs==None:
			l0_isogs=isogenies_2(E)
			E_to_l0_isogs[E]=l0_isogs
		

		for isog in l0_isogs:
			if testing:
				print(f"{format_j_path(j_path+[isog.codomain().j_invariant()])}")
			l0_h_isogs.append(hom_comp.from_factors(path+[isog]))
		nb_isogs+=base


		
		
	def go_back(base, h, intermediate_stop=_sage_const_0 ):
		nonlocal step, path_base_l0_plus_1, E


		while path_base_l0_plus_1[step]==base-_sage_const_1  and step>=intermediate_stop:
			path_base_l0_plus_1[step]=_sage_const_0 
			step-=_sage_const_1 
			E=path[step].domain()
			if args.verbose:
				print(f"At step {step} and curve {E.j_invariant()}")

		path_base_l0_plus_1[step]+=_sage_const_1 

		if step==intermediate_stop-_sage_const_1 :
			return True
		else:
			return False

	E=E0
	path_base_l0_plus_1=[_sage_const_1 ]*h
	l0_h_isogs=[]
	path=[]
	step=_sage_const_0 
	E_to_l0_isogs={}
	j_path=[]


	done=False
	nb_isogs=_sage_const_0 

	go_to_depth(h, path_base_l0_plus_1, init=True, testing=True)
	#Petit problème d'indice, calcul l0+1 fois la même chose avant de changer
	while not done:
		go_to_depth(h, path_base_l0_plus_1,			init=False, testing=True)

		add_isogs(l0+_sage_const_1 , path, testing=True)
		if args.verbose:
			print(f"Computed {nb_isogs} isogenies\n")
			print(f"Following path: {path_base_l0_plus_1}\n")
			print(f"Computed {l0+_sage_const_1 } leaves: {l0_h_isogs[nb_isogs-_sage_const_3 :nb_isogs]}\n")
		
		done=go_back(l0+_sage_const_1 , h,			intermediate_stop)


	return l0_h_isogs


	
	

#USE BRUTE_FORCE_ORIENTATION_SMALL instead
def brute_force_orientation(E0, primes, alpha1, alpha2, trace, norm):
	G=cartesian_product([[_sage_const_0 ,_sage_const_1 ,_sage_const_2 ]]*h)
	for seq in G:
		isog_factors=[]
		E_j=E0
	#should find a way to not recompute every l_0 isogs each time
		i=_sage_const_0 
		for j in seq:
			isogs=isogenies_2(E_j)
			isog_factors.append(isogs[j])
			E_j=isogs[j].codomain()
		#only compute kernel
		Phi=hom_comp.from_factors(isog_factors, E0)	
		print(Phi)
		print(Phi.degree())

	
w_0=compute_initial_orientation_l(E0, primes, alpha1, alpha2, trace, norm)

print(f"Computed initial orientation: {w_0} and factored degree {factor(w_0.degree())}")


"""isogs=brute_force_orientation_small(E0, primes, alpha1,alpha2,\
				trace, norm, l0, h,\
				intermediate_stop=Integer(args.samples))
"""

"""if args.verbose:
	print(f"Computed {len(isogs)} isogenies from {E0} of degree {l0}**{h}={l0**h}\n")
"""

F.close()

