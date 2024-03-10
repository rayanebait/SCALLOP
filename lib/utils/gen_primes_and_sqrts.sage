from sage.rings.finite_rings.integer_mod import square_root_mod_prime

K.<i>=NumberField(x^2+1)
nb_primes=11
primes=[]
p=5
for i in range(nb_primes):
	primes.append(p)
	p=p.next_prime()
	while p%4!=1:
		p=p.next_prime()
print(primes, len(primes))

sqrts=[square_root_mod_prime(Mod(-1,i), p=i) for i in primes]
print(sqrts, len(sqrts))
print([Mod(sqrts[i]**2, primes[i]) for i in range(nb_primes)])

filename="../../txt/sqrts_"+str(nb_primes)+"_primes.md"

f=open(filename, "w")

f.write(f"{nb_primes}\n")
for i in range(nb_primes):
	f.write(f"{sqrts[i]}\n")

f.close()
