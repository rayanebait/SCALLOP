
def next(L, primes, nb_primes, prev_seq, next_seq):
	for i in range(1, nb_primes):
		j=prev_seq[i]
		k=next_seq[i]
		if j!=k:
			L=L*(primes[i]**(4*k))
	return L

def gen_candidates(K, primes):
	nb_primes=len(primes)

	L_2=primes[0]
	L_1=prod(primes[1:len(primes)])**2
	L=L_1*L_2

	#could try not generating it each time
	G=cartesian_product([[-1,1] for i in range(nb_primes)])

	name="candidate_conductors"
	name=name+str(nb_primes)+".md"
	print(name)

	F=open(name, 'w', encoding="utf-8")
	F.write(f'{nb_primes}')
	F.write('\n')

	prev_seq=G[0]
	nb_candidates=0
	for next_seq in G:
		L=next(L,primes,nb_primes,prev_seq,next_seq)
		prev_seq=next_seq
		f_im=Integer(L.gens_reduced()[0].imag())
		f_re=Integer(L.gens_reduced()[0].real())
		
		pim=(abs(f_im)).is_pseudoprime()
		pre=(abs(f_re)).is_pseudoprime()
		if pim or pre:
			if pim:
				F.write(f'{f_im} {1}')
				F.write('\n')
				nb_candidates+=1
			else:
				F.write(f'{f_im} {0}')
				F.write('\n')

			if pre:
				F.write(f'{f_re} {1}')
				F.write('\n')
				nb_candidates+=1
			else:
				F.write(f'{f_re} {0}')
				F.write('\n')
	F.close()

def add_split_primes_below_B(primes, last_prime):
	last_prime=last_prime.next_prime()
	while last_prime%4!=1:
		last_prime=last_prime.next_prime()

	primes.append(last_prime)
	return (primes, last_prime)

def add_prime_above(p, last_prime):
	p.append(K.primes_above(last_prime)[0])
	return p

initial_B=90
B=100

K.<i>=NumberField(x^2+1)
last_prime=5
split_primes=[5]
p_above=[]

while last_prime<initial_B:
	(split_primes, last_prime)=add_split_primes_below_B(split_primes, last_prime)
	p_above=add_prime_above(p_above, last_prime)
print(last_prime)

while last_prime<B:
	print(f"next: {last_prime}\n")
	gen_candidates(K, p_above)
	print(p_above)
	print(split_primes)
	(split_primes, last_prime)=add_split_primes_below_B(split_primes, last_prime)
	p_above=add_prime_above(p_above, last_prime)
