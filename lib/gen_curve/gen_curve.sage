
def FullRepresentInteger(M, p, trials=0, s=0):
	m_=floor(float(sqrt(4*M/p)))
	z_=randrange(-m_, m_)
	s+=1
	seed(s)

	m__=floor(float( sqrt(4*M/p-z_**2) ))
	t_=randrange(-m__, m__)
	s+=1
	seed(s)

	M_=4*M-p*((z_)**2+(t_)**2)

	#Test if M_=1 mod 4
	is_sum=((M_&3)==1)
	if not is_sum:
		if trials>1000:
			return
		return FullRepresentInteger(M,p,trials+1, s)

	x_,y_=two_squares(M_)

	if (x_-t_)&1!=0 or (z_-y_)&1!=0:
		if trials>1000:
			return
		return FullRepresentInteger(M,p,trials+1, s)
	return (x_, y_, z_, t_)



L1=prod([p for p in range(5, 100) if p&3==1 and p in Primes()])**2
L2=5
L=L1*L2

p=1010000001983298392831832983928328938293829328932
f=3014688773870022715669219

h=floor(float(log(p/f)/log(97)))+1

M=(97**h)*f

print(float(sqrt(4*M/p)))

m_=floor(float(sqrt(4*M/p)))
z_=ZZ.random_element(-m_,m_)

m__=floor(float( sqrt(4*M/p-z_**2) ))
t_=ZZ.random_element(-m__,m__)

M_=4*M-p*((z_)**2+(t_)**2)

print(M_, z_, t_, m_, m__)

is_sum=(((M_)&3))
print(is_sum)
