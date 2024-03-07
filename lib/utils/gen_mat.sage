n=2
m=4

A = random_matrix(ZZ, n, m, x=-300, y=300)

A= Matrix(ZZ, 2,2, [[5, -8], [0, 1]])
print(A)

f=open("../txt/mat.md", 'w', encoding="utf-8")

f.write(f'{n} ')
f.write(f'{m} ')
for i in A.list():
	f.write(f' {i}')

f.close()
