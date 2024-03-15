from sage.rings.finite_rings.integer_mod import square_root_mod_prime
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-n', '--nbprimes', default='14')

args = parser.parse_args()

nb_primes=Integer(args.nbprimes)

K.<i>=NumberField(x^2+1)
primes=[]
p=5
i=0
while i<nb_primes:
	if p&3==1:
		primes.append(p)
		i+=1
	p=p.next_prime()

if args.verbose:
	print(f"Computing roots of -1 modulo the first {nb_primes}-primes\n")

sqrts=[square_root_mod_prime(Mod(-1,i), p=i) for i in primes]
if args.verbose:
	for i in range(nb_primes):
		print(f"{sqrts[i]}**2=-1 mod {primes[i]}\n")


filename="../../txt/sqrts_"+str(nb_primes)+"_primes.md"

f=open(filename, "w")

if args.verbose:
	print(f"Writing to file {filename}\n")

f.write(f"{nb_primes}\n")
for i in range(nb_primes):
	f.write(f"{sqrts[i]}\n")

f.close()
