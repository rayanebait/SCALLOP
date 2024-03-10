

# This file was *autogenerated* from the file gen_primes_and_sqrts.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_10 = Integer(10); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4)
from sage.rings.finite_rings.integer_mod import square_root_mod_prime

K = NumberField(x**_sage_const_2 +_sage_const_1 , names=('i',)); (i,) = K._first_ngens(1)
nb_primes=_sage_const_10 
primes=[]
p=_sage_const_5 
for i in range(nb_primes):
	primes.append(p)
	p=p.next_prime()
	while p%_sage_const_4 !=_sage_const_1 :
		p=p.next_prime()
print(primes, len(primes))

sqrts=[square_root_mod_prime(Mod(-_sage_const_1 ,i), p=i) for i in primes]
print(sqrts, len(sqrts))
print([Mod(sqrts[i]**_sage_const_2 , primes[i]) for i in range(nb_primes)])

filename="../../txt/primes_and_sqrts.md"
f=open(filename, "w")

f.write(f"{nb_primes}\n")
for i in range(nb_primes):
	f.write(f"{sqrts[i]}\n")

f.close()
