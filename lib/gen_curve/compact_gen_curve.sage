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

nb_primes=Integer(n)

if args.verbose:
	print(f"Generating oriented curve for {nb_primes} primes parameter\n")

endo_file="../../txt/endo_"+str(nb_primes)+"_primes.md"
prime_file="../../txt/prime_"+str(nb_primes)+"_primes.md"
cond_file="../../txt/conductor_"+str(nb_primes)+"_primes.md"
alpha_file="../../txt/w_"+str(nb_primes)+"_primes.md"

F=open(endo_file, "r")
G=open(prime_file, "r")
H=open(cond_file, "r")
O=open(alpha_file, "r")

p=Integer(G.readline())
K.<J>=GF((p,2), modulus=[1,0,1])
print(K)

conductor=Integer(O.readline())

alpha2=conductor
alpha1=Integer(O.readline())

norm=alpha1**2+alpha2**2
trace=2*alpha1

seq=[]
for i in range(nb_primes):
	seq.append(Integer(O.readline()))

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


def orientation_from_factors_slow(E_l, s_l):
	print(f"Generation of orientation for {E_l.j_invariant()}\n")

	(P_L1, P_L1_bar, P_L2)=s_l

	phi_L1=hom_comp(E_l, P_L1)
	phi_L1_bar=hom_comp(E_l, P_L1_bar).dual()

	E_L2=phi_L1.codomain()
	E_L2_bar=phi_L1_bar.domain()

	print(f"First isogeny of degree {factor(phi_L1.degree())}:\n\t{E_l.j_invariant()}--->{E_L2.j_invariant()}")

	phi_L2=E_L2.isogeny(phi_L1(P_L2))

	E_L2_bis=phi_L2.codomain()

	print(f"Second isogeny of degree {factor(phi_L2.degree())}:\n\t{E_L2.j_invariant()}--->{E_L2_bis.j_invariant()}")

	print(f"Third isogeny of degree {factor(phi_L1_bar.degree())}:\n\t{E_L2_bar.j_invariant()}--->{E_l.j_invariant()}\n")

	try:
		iso=E_L2_bis.isomorphism_to(E_L2_bar)
		w=phi_L1_bar*iso*phi_L2*phi_L1
	except ValueError:
		return None
	return w

def orientation_from_factors_try(E_l, s_l):
	print(f"Generation of orientation for {E_l.j_invariant()}\n")

	(P_L1, P_L1_bar, P_L2)=s_l

	phi_L1=hom_comp(E_l, P_L1)
	E_L2=phi_L1.codomain()
	
	phi_L2=E_L2.isogeny(phi_L1(P_L2))
	E_L2_bar=phi_L2.codomain()

	phi_L1_bar=hom_comp(E_l, P_L1_bar)
	E_L2_bis=phi_L1_bar.codomain()

	phi_L1_bar=hom_comp(E_L2_bis, phi_L1_bar(P_L1))


	print(f"First isogeny of degree {factor(phi_L1.degree())}:\n\t{E_l.j_invariant()}--->{E_L2.j_invariant()}")



	print(f"Second isogeny of degree {factor(phi_L2.degree())}:\n\t{E_L2.j_invariant()}--->{E_L2_bar.j_invariant()}")

	print(f"Third isogeny of degree {factor(phi_L1_bar.degree())}:\n\t{E_L2_bis.j_invariant()}--->{E_l.j_invariant()}\n")

	try:
		iso=E_L2_bar.isomorphism_to(E_L2_bis)
		iso_=phi_L1_bar.codomain().isomorphism_to(E_l)
		w=iso_*phi_L1_bar*iso*phi_L2*phi_L1
	except:
		return None
	return w
"""
	Compute the kernel of w_0 simply by getting kernels for 
	iota+/-sqrt_l mod l for each l. 

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
			seq,\
			keep_factors=False, testing=False):
	O=E0(0,0)-E0(0,0)

	gen_of_l_part=[]
	gen_of_l_part_bar=[]

	#Compute a generator of the kernel of each E[I_l]
	if args.verbose:
		print(f"Computing orientation w_0=alpha1+i*alpha2 with \n\
			\talpha1={alpha1}\n\
			\talpha2={alpha2}\n\
			And norm: {factor(norm)}")

	K1.<i>=NumberField(x^2+1)

	for (l, pmone) in zip(primes, seq):
		ind=0
		if pmone==-1:
			"""
			The ideals above l are of the form (l, i+sqrt_l), when 
			the conductor was chosen, it was chosen as a product of 
			primes above l's, pmone=-1 means second ideal above l
			"""
			ind=1
		else:
			ind=0

		gen=K1.primes_above(l)[ind].gens_two()[1]
		print(f"Generator for {l}: {gen}\n")
		sqrt_l=gen.real()
		if args.verbose:
			print(f"\nComputing generator for E[I] for I a prime above prime :{l}\n")
	
		P_l, Q_l=E0.torsion_basis(l)
		if args.verbose:
			print(f"Generators for the l-torsion:\n\t{P_l}, {Q_l}\n\n")

		if args.verbose:
			print(f"Pushing {P_l} and {Q_l} through iota+/-sqrt_l\n")


		"""It happens that w_0 has the same eigenspaces as iota"""
		(P_l_bar, Q_l_bar)=(sqrt_l*P_l+pmone*iota(P_l), sqrt_l*Q_l+pmone*iota(Q_l))
		(P_l, Q_l)=(sqrt_l*P_l+(-pmone)*iota(P_l), sqrt_l*Q_l+(-pmone)*iota(Q_l))

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

	if seq[0]==-1:
		ind=1
	else:
		ind=0

	P_L2=gen_of_l_part[ind]

	s_0 = (P_L1, P_L1_bar, P_L2)

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
	path_base_l0_plus_1=[0]*h
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
			trace, norm,\
			seq,\
			testing=testing)
	print(f"Computed initial orientation: {w_0}\n")
	exit(0)
else:
	s_0=compute_initial_orientation_l(\
			E0, primes,\
			alpha1, alpha2,\
			trace, norm, seq,\
			keep_factors=False)

if nb_primes==0:
	exit(0)
else:
	if args.verbose:
		print(f"Pushing initial orientation through endomorphism {endo[0]}+i*{endo[1]} of norm {factor(endo[0]**2+endo[1]**2)}\n")

	s=[endo[0]*P+endo[1]*(iota(P)) for P in s_0]
	
	

isogs=gen_isogs_small(E0, primes, alpha1, alpha2,\
				trace, norm, l0, h,\
				intermediate_stop=Integer(args.samples))


"""
We know tr(w) and 0 are eigenvalues of w mod l so that 
we just need to check if P and P_bar vanish 
through w*(w-tr)?
"""
def CheckTrace(Ei, w, tr):
	w_tr=abs(w.trace())
	print(f"Endomorphism w of trace {w_tr}, searched trace {tr}\n")
	if w_tr!=abs(tr):
		print(f"Checking trace: {w.trace()==tr}")
		return 0
	return 1
		

w=iota
E_l=E0
for phi_i in isogs:
		
	Ei=phi_i.codomain()
	s_Ei=[phi_i(P) for P in s]

	w_Ei=orientation_from_factors_try(Ei, s_Ei)
	if w_Ei==None:
		continue

	if w_Ei.degree()!=norm:
		continue

	if CheckTrace(Ei, w_Ei, trace)==1:
		w=w_Ei
		E_l=phi_i.codomain()
		break
	

if w==iota or E_l==E0:
	print(f"Curve generation failed, didn't find the matching l_0^h isogeny\n")
else:
	print(f"Successfully generated orientation for {E_l.j_invariant()}:\n\t {w}\n")
	print(f"With trace: {w.trace()}")

F.close()
