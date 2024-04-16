from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-in','--initial_nb_primes', default='2')
parser.add_argument('-mn','--maximal_nb_primes', default='3')
parser.add_argument('-v','--verbose', action='store_true')

args=parser.parse_args()

def next(L, p_above,p_above_conj, nb_primes, prev_seq, next_seq):
	for i in range(nb_primes-1):
		j=prev_seq[i]
		k=next_seq[i]
		if j!=k:
			L=L*(p_above[i+1]**(2*k))
			L=L*(p_above_conj[i+1]**(2*(-k)))
	return L

def gen_candidates(K, p_above, p_above_conj, f_abs, a_abs, stop=1000000, G_index=0):
	nb_primes=len(p_above)
	print(nb_primes)

	L_2=p_above[0]
	L_1=prod(p_above_conj[1:len(p_above)])**(2)
	L=L_1*L_2

	#could try not generating it each time
	G=cartesian_product([[-1,1] for i in range(nb_primes-1)])


	prev_seq=G[0]
	i=0
	for next_seq in G:
		if i<G_index:
			i+=1
			continue
		if stop==0:
			break

		L=next(L,p_above, p_above_conj,nb_primes,prev_seq,next_seq)
		prev_seq=next_seq

		f_im=Integer(L.gens_reduced()[0].imag())
		f_re=Integer(L.gens_reduced()[0].real())

		if args.verbose:
			print(f"Real part of w: {f_re}")
			print(f"Imaginary part of w: {f_im}")

		if args.verbose:
			print(f"Corresponding sequence: {[1]+list(prev_seq)}")

		if abs(f_re) in [f_abs, a_abs] and abs(f_re) in [f_abs, a_abs]:
			which=0
			if abs(f_re)==a_abs:
				which=1
			return (f_re, f_im, [1]+list(prev_seq), which)
		i+=1


def add_prime_above(p_above, p_above_conj, p):
	primes_above=K.primes_above(p)
	p_above.append(primes_above[0])
	p_above_conj.append(primes_above[1])
	return (p_above, p_above_conj)


initial_nb_primes=Integer(args.initial_nb_primes)
max_nb_primes=Integer(args.maximal_nb_primes)

K.<i>=NumberField(x^2+1)
p=5
split_primes=[]
p_above=[]
p_above_conj=[]
i=0

while i<initial_nb_primes:
	if p&3==1:
		split_primes.append(p)
		(p_above, p_above_conj)=add_prime_above(p_above, p_above_conj, p)
		i+=1
	p=p.next_prime()

nb_primes=max_nb_primes

while nb_primes<=max_nb_primes:
	while p&3!=1:
		p=p.next_prime()

	split_primes.append(p)
	(p_above, p_above_conj)=add_prime_above(p_above, p_above_conj, p)

	if args.verbose:
		print(f"Searching w for primes {split_primes}")

	name1="../../txt/alpha_"+str(nb_primes)+"_primes.md"
	name2="../../txt/w_"+str(nb_primes)+"_primes.md"
	if args.verbose:
		print(f"Reading alpha in {name1}")

	F=open(name1, 'r', encoding="utf-8")
	f_abs=Integer(F.readline())
	a_abs=Integer(F.readline())
	F.close()
	if args.verbose:
		print(f"alpha={a_abs}+i{f_abs}")

	res=gen_candidates(K, p_above, p_above_conj, f_abs, a_abs, stop=100000)
	if args.verbose:
		print(f"Found corresponding w {res[0]}+i{res[1]}\n corresponding to sequence\n\t{res[2]}\n")

	if args.verbose:
		print(f"Writing w in {name2}")

	G=open(name2, 'w', encoding="utf-8")
	G.write(str(res[res[3]]))
	G.write('\n')
	G.write(str(res[abs(res[3]-1)]))
	G.write('\n')
	seq=res[2]
	for pmone in seq:
		G.write(str(pmone))
		G.write('\n')
	G.close()


	p=p.next_prime()
	nb_primes+=1


