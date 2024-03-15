conductor=40
endo_file="../../txt/endo_"+str(conductor)+"_bits.md"
prime_file="../../txt/prime_"+str(conductor)+"_bits.md"
cond_file="../../txt/conductor_"+str(conductor)+"_bits.md"

F=open(endo_file, "r")
G=open(prime_file, "r")
H=open(cond_file, "r")

K=GF(p**2)

E=EllipticCurve()
