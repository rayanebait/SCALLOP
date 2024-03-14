primes=[]
p=5
i=0
while i<10:
	if p&3==1:
		primes.append(p)
	p=p.next_prime()
L1=floor(prod(primes)/5)
L2=5

f=3014688773870022715669219

(le L3 est peut etre la pour avoir une sécurité plus grosse)
p=L-1
c=1
is_prime=False
while True:
	p=c*L-1
	is_prime=p.is_pseudoprime()
	if p.is_pseudoprime():
		break
	p=p+2
	if p.is_pseudoprime():
		break
	c+=1

