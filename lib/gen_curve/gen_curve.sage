def FullRepresentInteger(M, p, trials=0):
	m_=floor(float(sqrt(4*M/p)))
	z_=ZZ.random_element(-m_, m_)

	m__=floor(float( sqrt( (4*M/p) - z_**2) ))
	t_=ZZ.random_element(-m__, m__)

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


conductor=80
nb_primes=17

cond_file="../../txt/conductor_"+str(conductor)+"_bits.md"
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

prime_file="../../txt/prime_"+str(conductor)+"_bits.md"
G=open(prime_file, "r")
p=Integer(G.readline())
G.close()



print(p,f)
l_0=2
while l_0<100:
	h=floor(float(log(p/f)/log(l_0)))+1
	M=(l_0**h)*f

	endo=FullRepresentInteger(M,p)

	if endo==None:
		l_0=l_0.next_prime()
		while l_0&3==1:
			l_0=l_0.next_prime()
	else:
		break

print(endo, l_0, h)

endo_file="../../txt/endo_"+str(conductor)+"_bits.md"
K=open(endo_file,"w")
for coeff in endo:
	K.write(f"{coeff}\n")
K.write(f"{l_0}\n{h}\n")
K.close()
