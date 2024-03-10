filename="../../txt/conductor_80_bits.md"

f=open(filename, "r")
f.readline()

conductor=Integer(f.readline())

f.close()

filename="../ideals/candidate_conductors16.md"
g=open(filename, "r")
line=0

n=abs(Integer(g.readline().split(" ")[0]))

while n!=conductor:
	n=abs(Integer(g.readline().split(" ")[0]))
	line+=1

print(line)
print(conductor,n)
