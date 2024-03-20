from argparse import ArgumentParser
parser=ArgumentParser()

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-n', '--nbprimes', default='14')

args=parser.parse_args()

proof.arithmetic(False)

def FullRepresentInteger(M, p, trials=0):
	m_=floor(float(sqrt(4*M/p)))
	z_=ZZ.random_element(-m_, m_+1)

	m__=floor(float( sqrt( (4*M/p) - z_**2) ))
	t_=ZZ.random_element(-m__, m__+1)

	M_=4*M-p*((z_)**2+(t_)**2)

	try:
		x_,y_=two_squares(M_)
	except ValueError:
		if trials>100:
			return
		return FullRepresentInteger(M,p,trials+1)

	if (x_-t_)&1!=0 or (z_-y_)&1!=0:
		if trials>100:
			return
		return FullRepresentInteger(M,p,trials+1)
	return (x_, y_, z_, t_)

nb_primes=Integer(args.nbprimes)

cond_file="../../txt/conductor_"+str(nb_primes)+"_primes.md"
if args.verbose:
	print(f"Reading conductor from {cond_file}\n")

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

L1=prod(primes[1:])**2
L2=5
L=L1*L2

prime_file="../../txt/prime_"+str(nb_primes)+"_primes.md"
if args.verbose:
	print(f"Reading characteristic from {prime_file}\n")

G=open(prime_file, "r")
p=Integer(G.readline())
G.close()


if args.verbose:
	print(f"Startng endomorphism generation with FullRepresentInteger\n")
l_0=5
while l_0<100:
	if args.verbose:
		print(f"Trying with l_0={l_0}")

	h=floor(float(log(p/f)/log(l_0)))+1
	M=(l_0**h)*f

	endo=FullRepresentInteger(M,p)

	if endo==None:
		l_0=l_0.next_prime()
		while l_0&3==1:
			l_0=l_0.next_prime()
	else:
		break

if args.verbose:
	print(f"Found endomorphism {endo}/2 of norm M=(l_0**h)*f={l_0}**{h}*{f}\n")

endo_file="../../txt/endo_"+str(nb_primes)+"_primes.md"
if args.verbose:
	print(f"Writing endomorphism to file {endo_file}\n")

K=open(endo_file,"w")
for coeff in endo:
	K.write(f"{coeff}\n")
K.write(f"{l_0}\n{h}\n")
K.close()
