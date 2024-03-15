nb_primes=17
conductor=80

primes=[]
p=5
i=0
while i<nb_primes:
	if p&3==1:
		primes.append(p)
		i+=1
	p=p.next_prime()
print(primes)

L1=prod(primes[1:])**2
L2=5
L=L1*L2

f=3014688773870022715669219
#f_squared=f**2

#(le L3 est peut etre la pour avoir une sécurité plus grosse)

print(L1, L2)
p=L-1
c=1
pmone=1

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
	print(c)
print(p in Primes())
print(p, c, pmone)

filename="../../txt/prime_"+str(conductor)+"_bits.md"
G=open(filename,"w")

G.write(f"{p}\n")
G.write(f"{c}\n{pmone}")
G.close()
