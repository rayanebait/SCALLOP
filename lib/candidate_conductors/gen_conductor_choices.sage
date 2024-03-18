def next(L, p_above,p_above_conj, nb_primes, prev_seq, next_seq):
	#print(prev_seq, next_seq)
	for i in range(nb_primes-1):
		j=prev_seq[i]
		k=next_seq[i]
		if j!=k:
			L=L*(p_above[i+1]**(2*k))
			L=L*(p_above_conj[i+1]**(2*(-k)))
	return L

def gen_candidates(K, p_above, p_above_conj, stop=1000000, G_index=0):
	nb_primes=len(p_above)
	print(nb_primes)

	L_2=p_above[0]
	L_1=prod(p_above_conj[1:len(p_above)])**(2)
	L=L_1*L_2

	#could try not generating it each time
	G=cartesian_product([[-1,1] for i in range(nb_primes-1)])

	name="candidate_conductors"+str(nb_primes)+".md"
	print(name)

	F=open(name, 'w', encoding="utf-8")
	F.write(f'{nb_primes}\n')

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
		
		pim=(abs(f_im)).is_pseudoprime()
		pre=(abs(f_re)).is_pseudoprime()
		if pim or pre:
			if pim:
				F.write(f'{f_im} {1}')
				F.write('\n')
				stop-=1
			else:
				F.write(f'{f_im} {0}')
				F.write('\n')

			if pre:
				F.write(f'{f_re} {1}')
				F.write('\n')
				stop-=1
			else:
				F.write(f'{f_re} {0}')
				F.write('\n')
		i+=1
	F.write(f"{i}")
	F.close()


def add_prime_above(p_above, p_above_conj, p):
	primes_above=K.primes_above(p)
	p_above.append(primes_above[0])
	p_above_conj.append(primes_above[1])
	return (p_above, p_above_conj)


initial_nb_primes=13
max_nb_primes=14

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

nb_primes=initial_nb_primes

while nb_primes<max_nb_primes:
	while p&3!=1:
		p=p.next_prime()

	split_primes.append(p)
	print(f"\n next: {p}\n")
	(p_above, p_above_conj)=add_prime_above(p_above, p_above_conj, p)

	gen_candidates(K, p_above, p_above_conj, stop=100000)
	print(split_primes)

	p=p.next_prime()
	nb_primes+=1
