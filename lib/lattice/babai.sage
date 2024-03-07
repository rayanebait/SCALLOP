#A a pour colonnes la base, A_ est la base réduite
#t est le vecteur dont on cherche le vecteur proche

def babai(A,A_, t):
	s=-t

	I=A.columns()
	I.reverse()

	J=A_.columns()
	J.reverse()

	vecs=zip(I,J)
	c=1
	for (v,v_) in vecs:
		c=floor((s*v_)/(v_*v_))
		s=s-c*v
	#In the gram_schmidt lattice
	return s

e=10**9
size=3*floor(float(log(e)))

print("nombre de vecteurs dans la base:", size)

M=MatrixSpace(ZZ,size)
#A=diagonal_matrix(ZZ, size, vector(ZZ, [ZZ.random_element(1, 1000) for i in 
#range(size)]))
A = random_matrix(ZZ, size, size, x=-300, y=300)
#A=M.random_element()


while A.rank() != size:
#ensures V has rank size
	A=M.random_element()
	print(A)



t=vector(ZZ, [e]+(size-1)*[0])
#t=A.columns()[0]*10000

print("réduction d'un vecteur de taille:", float(log(e)))

A=A.BKZ(block_size=floor(sqrt(size)), proof=True)
V=span(A.columns(), ZZ)

G,_=A.gram_schmidt()
#print(G)

t_2=babai(A,G,t)
print(float(log(max(t_2))), (t_2+t) in V)
t_3=babai(A,G,t_2)
print(float(log(max(t_3))), (t_3+t_2) in V)
t_4=babai(A,G,t_3)
print(float(log(max(t+t_4))), (t_4+t_3) in V)
