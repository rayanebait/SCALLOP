

# This file was *autogenerated* from the file gen_conductor_choices.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_250 = Integer(250); _sage_const_260 = Integer(260); _sage_const_5 = Integer(5)
def next(L, primes, nb_primes, prev_seq, next_seq):
	for i in range(_sage_const_1 , nb_primes):
		j=prev_seq[i]
		k=next_seq[i]
		if j!=k:
			L=L*(primes[i]**(_sage_const_4 *k))
	return L

def gen_candidates(K, primes):
	nb_primes=len(primes)

	L_2=primes[_sage_const_0 ]
	L_1=prod(primes[_sage_const_1 :len(primes)])**_sage_const_2 
	L=L_1*L_2

	#could try not generating it each time
	G=cartesian_product([[-_sage_const_1 ,_sage_const_1 ] for i in range(nb_primes)])

	name="candidate_conductors"
	name=name+str(nb_primes)+".md"
	print(name)

	F=open(name, 'w', encoding="utf-8")
	F.write(f'{nb_primes}')
	F.write('\n')

	prev_seq=G[_sage_const_0 ]
	nb_candidates=_sage_const_0 
	for next_seq in G:
		L=next(L,primes,nb_primes,prev_seq,next_seq)
		prev_seq=next_seq
		f_im=Integer(L.gens_reduced()[_sage_const_0 ].imag())
		f_re=Integer(L.gens_reduced()[_sage_const_0 ].real())
		
		pim=(abs(f_im)).is_pseudoprime()
		pre=(abs(f_re)).is_pseudoprime()
		if pim or pre:
			if pim:
				F.write(f'{f_im} {_sage_const_1 }')
				F.write('\n')
				nb_candidates+=_sage_const_1 
			else:
				F.write(f'{f_im} {_sage_const_0 }')
				F.write('\n')

			if pre:
				F.write(f'{f_re} {_sage_const_1 }')
				F.write('\n')
				nb_candidates+=_sage_const_1 
			else:
				F.write(f'{f_re} {_sage_const_0 }')
				F.write('\n')
	F.close()

def add_split_primes_below_B(primes, last_prime):
	last_prime=last_prime.next_prime()
	while last_prime%_sage_const_4 !=_sage_const_1 :
		last_prime=last_prime.next_prime()

	primes.append(last_prime)
	return (primes, last_prime)

def add_prime_above(p, last_prime):
	p.append(K.primes_above(last_prime)[_sage_const_0 ])
	return p

initial_B=_sage_const_250 
B=_sage_const_260 

K = NumberField(x**_sage_const_2 +_sage_const_1 , names=('i',)); (i,) = K._first_ngens(1)
last_prime=_sage_const_5 
split_primes=[_sage_const_5 ]
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

