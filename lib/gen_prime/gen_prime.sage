from argparse import ArgumentParser
parser=ArgumentParser()

parser.add_argument('-v','--verbose',action='store_true')
parser.add_argument('-n','--nbprimes', default='14')

args=parser.parse_args()
nb_primes=Integer(args.nbprimes)

cond_file="../../txt/conductor_"+str(nb_primes)+"_primes.md"
F=open(cond_file, "r")
F.readline()

f=Integer(F.readline())

F.close()
primes=[]
p=5
i=0
while i<nb_primes:
	if p&3==1:
		primes.append(p)
		i+=1
	p=p.next_prime()


L1=prod(primes[1:])
L2=5
L=L1*L2

p=L-1
c=1
pmone=1

if args.verbose:
	print(f"Generating characteristic for conductor {f} with norm {L} and primes {primes}\n")
while True:
	p=c*L-1
	if p.is_pseudoprime() and p&3==3:
		pmone=-1
		break
	p=p+2
	if p.is_pseudoprime() and p&3==3:
		pmone=1
		break
	c+=1

if args.verbose:
	if p in Primes():
		print(f"Found prime {p} with c={c} and p=cL+({pmone})\n")
	else:
		print("Pseudo-primality test failed\n")
		raise ValueError

filename="../../txt/prime_"+str(nb_primes)+"_primes.md"
if args.verbose:
	print(f"Storing prime {p} in file {filename}")

G=open(filename,"w")
G.write(f"{p}\n{c}\n{pmone}")
G.close()
