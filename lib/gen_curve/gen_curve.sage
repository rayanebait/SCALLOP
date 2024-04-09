from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from sage.schemes.elliptic_curves.weierstrass_morphism import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite as hom_comp
from sage.schemes.elliptic_curves.isogeny_small_degree import isogenies_2
from argparse import ArgumentParser
import time 
from pprint import pp
from os import get_terminal_size

proof.arithmetic(False)

parser = ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-n', '--nbprimes', default='3')
parser.add_argument('-s', '--samples', default='0')
parser.add_argument('-t', '--testing', action='store_true')
args = parser.parse_args()


n=args.nbprimes
if n!='3' and n!='14' and n!='19' and n!='20' and n!='26':
	raise SystemExit(f"Currently supports 3, 14, 19, 20 and 26 primes parameters, not {n}\n")

nb_lines_terminal=get_terminal_size().lines
CLEAR='\r'+'\033[<'+str(nb_lines_terminal)+'>A'

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

def orientation_from_factors(E_l, l_parts_of_s_l):
	E_i, E_i_bar=(E_l,E_l)
	(gen_of_l_part, gen_of_l_part_bar)=l_parts_of_s_l

	j_path=[E_l.j_invariant()]
	j_bar_path=[E_l.j_invariant()]

	isogeny_chain=[]
	isogeny_chain_bar=[]
	curve_chain=[]
	curve_chain_bar=[]


	#compute the isogeny factors of phi_L1 and phi_L1^-1
	for (P,Pbar) in zip(gen_of_l_part[1:], gen_of_l_part_bar[1:]):
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
	
	P_L2=gen_of_l_part[0]
	P_L2_bar=gen_of_l_part_bar[0]
	for (phi,phi_bar) in zip(isogeny_chain, isogeny_chain_bar):
		P_L2=phi(P_L2)
		P_L2_bar=phi_bar(P_L2_bar)

	phi_L2=E_i.isogeny(P_L2)
	phi_L2_bar=E_i.isogeny(P_L2_bar)

	E_L2=E_i

	E_L2_bis=phi_L2.codomain()
	E_L2_bis_bis=phi_L2_bar.codomain()

	phi_L1=hom_comp.from_factors(isogeny_chain, E_l)
	phi_L1bar=(hom_comp.from_factors(isogeny_chain_bar, E_l)).dual()

	if args.verbose:
		path_str='--->'.join(str(j) for j in j_path[:len(j_path)-1])
		path_str+="--->"+str(j_path[-1])

		print(f"Took path:\n\t")
		print(f"{path_str} for phi_L1\n\n")

		path_str=str(E_L2.j_invariant())+"--->"+str(E_L2_bis.j_invariant())
		path_str_bar=str(E_L2.j_invariant())+"--->"+str(E_L2_bis_bis.j_invariant())

		print(f"Took path:\n\t")
		print(f"{path_str} for phi_L2\n\n")
		print(f"and path:\n\t")
		print(f"{path_str_bar} for phi_L2_bar\n\n")

		path_str='--->'.join(str(j) for j in j_bar_path[:len(j_bar_path)-1])
		path_str+="--->"+str(j_bar_path[-1])

		print(f"Took path:\n\t")
		print(f"{path_str} for phi_L1^-1\n\n")

	try:
		iso=E_L2_bis.isomorphism_to(E_i_bar)
		w=phi_L1bar*iso*phi_L2*phi_L1
	except TypeError:
		try:
			iso=E_L2_bis_bis.isomorphism_to(E_i_bar)
			w=phi_L1bar*iso*phi_L2*phi_L1
		except TypeError:
			raise RuntimeError(f"Orientation generation failed\n")

	return w

def orientation_from_factors_slow(E_l, s_l):
	print(f"Generation of orientation for {E_l.j_invariant()}\n")

	(P_L1, P_L1_bar, (P_L2, P_L2_bar))=s_l

	phi_L1=hom_comp(E_l, P_L1)
	phi_L1_bar=hom_comp(E_l, P_L1_bar).dual()

	E_L2=phi_L1.codomain()
	E_L2_bar=phi_L1_bar.domain()

	print(f"First isogeny of degree {factor(phi_L1.degree())}:\n\t{E_l.j_invariant()}--->{E_L2.j_invariant()}")

	phi_L2=E_L2.isogeny(phi_L1(P_L2))
	phi_L2_bar=E_L2.isogeny(phi_L1(P_L2_bar))

	E_L2_bis=phi_L2.codomain()
	E_L2_bis_bis=phi_L2_bar.codomain()

	print(f"Second isogeny of degree {factor(phi_L2.degree())} between:\n\t{E_L2.j_invariant()}--->{E_L2_bis.j_invariant()}")
	print(f"and (of degree {factor(phi_L2_bar.degree())}):\n\t{E_L2.j_invariant()}--->{E_L2_bis_bis.j_invariant()}")

	print(f"Third isogeny of degree {factor(phi_L1_bar.degree())}:\n\t{E_L2_bar.j_invariant()}--->{E_l.j_invariant()}\n")

	try:
		iso=E_L2_bis.isomorphism_to(E_L2_bar)
		w=phi_L1bar*iso*phi_L2*phi_L1
	except ValueError:
		None

	try:
		iso=E_L2_bis_bis.isomorphism_to(E_L2_bar)
		w=phi_L1bar*iso*phi_L2*phi_L1
	except ValueError:
		return None
	return w

"""
	Compute the kernel of w_0 simply by getting kernels for 
	w_0-Tr(w_0) mod l for each l. This works because (deg(w_0)=0
	mod l so that (w_0-Tr(w_0))*w_0=0). 
f

	To build the morphism choose an ordering on the l's, say
	(l_i)_i and push l_i through the chain:

		phi_l_(i-1)*...*phi_l_0

	to get phi_l_i. Do the same for l_i_bar without l_0 to
	complete the lollipop. (isomorphisms may need to be taken
	to compose phi_L1L2 and the dual of phi_L1^-1.
"""
def compute_initial_orientation_l(\
			E0, primes, alpha1,\
			alpha2,trace, norm,\
			keep_factors=False, testing=False):
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

	for l in primes:
		sqrt_l=Integer(square_root_mod_prime(Mod(-1, l),p=l))
		if args.verbose:
			print(f"\nComputing generator for E[I] for I a prime above prime :{l}\n")
	
		P_l, Q_l=E0.torsion_basis(l)
		if args.verbose:
			print(f"Generators for the l-torsion:\n\t{P_l}, {Q_l}\n\n")

		if args.verbose:
			print(f"Pushing {P_l} and {Q_l} through iota+/-sqrt_l\n")


		(P_l_bar, Q_l_bar)=(sqrt_l*P_l-iota(P_l), sqrt_l*Q_l-iota(Q_l))
		(P_l, Q_l)=(sqrt_l*P_l+iota(P_l), sqrt_l*Q_l+iota(Q_l))

		"""Now P_l, Q_l\in ker(iota-sqrt_l), they even generate it if nonzero (one of them must be)."""

		if P_l==O:
			gen_of_l_part.append(Q_l)
		else:
			gen_of_l_part.append(P_l)

		if P_l_bar==O:
			gen_of_l_part_bar.append(Q_l_bar)
		else:
			gen_of_l_part_bar.append(P_l_bar)

	P_L1=sum(gen_of_l_part[1:])
	P_L1_bar=sum(gen_of_l_part_bar[1:])
	gens_L2=(gen_of_l_part[0], gen_of_l_part_bar[0])

	print(f"Orders: {factor(P_L1.order())}, {factor(P_L1_bar.order())}, {gens_L2[0].order()}, {gens_L2[1].order()}\n")

	s_0 = (P_L1, P_L1_bar, gens_L2)

	if testing:
		return orientation_from_factors_slow(E0, s_0)
	else:
		if args.verbose:
			print(f"Found kernel representation {s_0} for initial orientation\n")
		if keep_factors:
			return (gen_of_l_part, gen_of_l_part_bar)
		else:
			return s_0

	
endo=[]
for i in range(4):
	endo.append(Integer(F.readline())//2)

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



def gen_isogs_small(E0, primes, alpha1,\
				alpha2, trace, norm,\
				l0, h,\
				intermediate_stop=0):

	if intermediate_stop>h :
		raise ValueError(f"-i should be < to h\n")
		return
	elif intermediate_stop<h-10:
		"""For at most l0^10 samples at once, should make the 
		function yield l0^10 isogeny at once, just to make it lighter
		"""
		if h<10:
			intermediate_stop=0
		else:
			intermediate_stop=h-10

	def format_j_path(j_path):
		return "--->".join(str(j) for j in j_path[:len(j_path)-1])\
			+"--->"+str(j_path[-1])

	def go_to_depth(base, h, path_base_l0_plus_1,\
			init=False, testing=False):

		nonlocal path, step,\
			E_to_l0_isogs, E_i,\
			l0_h_isogs, nb_isogs

		if testing:
			nonlocal j_path

		if step==h-1:
			if testing:
				if init:
					j_path.append(E_i.j_invariant())
				else:
					j_path[step]=E_i.j_invariant()
				print(f"{format_j_path(j_path)}\nwith base {base} path {path_base_l0_plus_1} and length {h}\n")


			l0_h_isogs.append(\
				hom_comp.from_factors(path)\
			)
			nb_isogs+=1

			return

		while step < h:
			if testing:
				print(f"At node {E_i.j_invariant()} at step {step}\n")
			"""
			Check if the codomain has already been encountered if yes,
			get the l0-isogenies. Then at each step, check if the 
			isogeny chosen is the dual of the last one, if yes jump it.
			(seems not doable since dual is defined up to isomorphism)
			"""
			l0_isogs=E_to_l0_isogs.get(E_i)
			if l0_isogs == None:
				l0_isogs=isogenies_2(E_i)
				E_to_l0_isogs[E_i]=l0_isogs
			isog=l0_isogs[path_base_l0_plus_1[step]]

			if init:
				path.append(isog)
				if testing:
					j_path.append(E_i.j_invariant())
			else: 
				path[step]=isog
				if testing:
					j_path[step]=E_i.j_invariant()
			
			E_i=path[step].codomain()
			step+=1
		step-=1
		if testing:
			if init:
				j_path.append(E_i.j_invariant())
			else:
				j_path[step]=E_i.j_invariant()
			print(f"{format_j_path(j_path)}\nwith base {base} path {path_base_l0_plus_1} and length {h}\n")


		l0_h_isogs.append(\
			hom_comp.from_factors(path)\
		)
		nb_isogs+=1

		return

	def add_isogs(base, path, testing=False):
		nonlocal E_i, E_to_l0_isogs,\
			nb_isogs, j_path,\
			h, step

		if testing:
			nonlocal j_path

		l0_isogs=E_to_l0_isogs.get(E_i)
		if l0_isogs==None:
			l0_isogs=isogenies_2(E_i)
			E_to_l0_isogs[E_i]=l0_isogs
		

		for isog in l0_isogs:
			if isog==path[step-1].dual():
				continue
			else:
				l0_h_isogs.append(\
					hom_comp.from_factors(path+[isog], strict=True)\
				)
			if testing:
				print(f"{format_j_path(j_path+[isog.codomain().j_invariant()])}\nwith base {base} path {path_base_l0_plus_1} and length {h}\n")
		nb_isogs+=base-1


		
		
	def go_back(base, h, intermediate_stop=0, testing=False):
		nonlocal step, path_base_l0_plus_1, E_i


		while path_base_l0_plus_1[step]==base-1 and step>=intermediate_stop:
			path_base_l0_plus_1[step]=0
			step-=1
			E_i=path[step].domain()
			if testing:
				print(f"At step {step} and curve {E_i.j_invariant()}")

		path_base_l0_plus_1[step]+=1

		if step==intermediate_stop-1:
			return True
		else:
			return False

	E_i=E0
	path_base_l0_plus_1=([1]*h)
	l0_h_isogs=[]
	path=[]
	step=0
	E_to_l0_isogs={}
	j_path=[]


	done=False
	nb_isogs=0

	go_to_depth(l0+1, h, path_base_l0_plus_1,\
			init=True, testing=False)

	while not done:
		go_to_depth(l0+1, h, path_base_l0_plus_1,\
			init=False, testing=False)
		
		done=go_back(l0+1, h,\
			intermediate_stop, testing=False)
	return l0_h_isogs
	

"""USE GEN_ISOGS_SMALL instead"""
def gen_isogs(E0, primes, alpha1, alpha2, trace, norm):
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

	
testing=args.testing
keep_factors=False
if testing:
	w_0=compute_initial_orientation_l(\
			E0, primes,\
			alpha1, alpha2,\
			trace, norm, testing)
	print(f"Computed initial orientation: {w_0}\n")
	exit(0)
else:
	if keep_factors:
		l_parts_of_s_0=compute_initial_orientation_l(\
				E0, primes,\
				alpha1, alpha2,\
				trace, norm, keep_factors=True)
	else:
		s_0=compute_initial_orientation_l(\
				E0, primes,\
				alpha1, alpha2,\
				trace, norm, keep_factors=False)

if nb_primes==0:
	exit(0)
else:
	if args.verbose:
		print(f"Pushing initial orientation through endomorphism {endo[0]}+i*{endo[1]} of norm {factor(endo[0]**2+endo[1]**2)}\n")

	if keep_factors:
		l_parts_of_s=([endo[0]*P_l+endo[1]* iota(P_l) \
				for P_l in l_parts_of_s_0[0]],\
			[endo[0]*P_l_bar+endo[1]*iota(P_l_bar) \
				for P_l_bar in l_parts_of_s_0[1]])
	else:
		s=(endo[0]*s_0[0]+endo[1]*iota(s_0[0]),\
			endo[0]*s_0[1]+endo[1]*iota(s_0[1]),\
			(endo[0]*P_L2+endo[1]*iota(P_L2) for P_L2 in s_0[2]))
	
	

isogs=gen_isogs_small(E0, primes, alpha1, alpha2,\
				trace, norm, l0, h,\
				intermediate_stop=Integer(args.samples))


"""
We know tr(w) and 0 are eigenvalues of w mod l so that 
we just need to check if P and P_bar vanish 
through w*(w-tr)?
"""
def CheckTrace(w, tr, l_parts_of_w):
	O=l_parts_of_w[0][0]-l_parts_of_w[0][0]

	for (P, P_bar) in zip(l_parts_of_w[0], l_parts_of_w[1]):
		(P, P_bar)=(w(w(P)-tr*P), w(w(P_bar)-tr*P_bar))
		if P!=O or P_bar!=O:
			return 0
	return 1
		
w=iota
E=E0


for phi_i in isogs:
	Ei=phi_i.codomain()
	if keep_factors:
		l_parts_of_s_Ei=([phi_i(P_l) for P_l in l_parts_of_s[0]],\
						[phi_i(P_l_bar) for P_l_bar in l_parts_of_s[1]])
		w1=orientation_from_factors(Ei, l_parts_of_s_Ei)
	else:
		s_Ei=(phi_i(s[0]),\
		   phi_i(s[1]),\
		   (phi_i(P) for P in s[2]))
		w1=orientation_from_factors_slow(Ei, s_Ei)
		if w1==None:
			break
	print(w1)
	continue
	"""
	Should 1728 be allowed ? 
	"""
	if Ei.j_invariant()==1728:
		continue

	"""Seems really slow, maybe just keep the factors"""
	phi_L1L2=hom_comp(Ei, s_Ei[0])
	phi_L1bar=hom_comp(Ei, s_Ei[1]).dual()

	"""Should check if the orientation can really be computed that way"""
	E=phi_L1L2.codomain()
	E_=phi_L1bar.domain()
	print("--->".join( [str(E.j_invariant()),str(E_.j_invariant())] ))
	try:
		iso=E.isomorphism_to(E_)
	except ValueError:
		continue

	w_Ei=phi_L1bar*iso*phi_L1L2
	print(factor(w_Ei.degree()))

	if w_Ei.degree()!=norm:
		continue
	l_parts_of_w_Ei=([phi_i(P_l) for P_l in l_parts_of_s[0]],\
			[phi_i(P_l_bar) for P_l_bar in l_parts_of_s[1]])

	if CheckTrace(w_Ei, trace, l_parts_of_w_Ei)==1:
		w=w_Ei
		E=phi_i.codomain()
		break
	

if w==iota or E==E0:
	print(f"Curve generation failed, didn't find the matching l_0^h isogeny\n")
else:
	print(f"Successfully generated orientation for {E.j_invariant()}:\n\t {w}\n")

F.close()
