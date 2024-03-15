from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-p', '--parameter', choices=['40', '80'], default='40')

args = parser.parse_args()

if args.parameter=="128":
	raise SystemExit("Aborting, this parameter size will take forever")
if args.parameter!="40" and args.parameter!="80":
	raise RuntimeError(f"Unsupported parameter {args.parameter}. Currently supports 40 and 80 bit parameters\n")

conductor=Integer(args.parameter)

if args.verbose:
	print(f"Generating oriented curve for {conductor}-bit security\n")

endo_file="../../txt/endo_"+str(conductor)+"_bits.md"
prime_file="../../txt/prime_"+str(conductor)+"_bits.md"
cond_file="../../txt/conductor_"+str(conductor)+"_bits.md"

F=open(endo_file, "r")
G=open(prime_file, "r")
H=open(cond_file, "r")

p=Integer(G.readline())

K.<i>=GF((p,2), modulus=[1,0,1])
print(K, i**2)


F.close()
G.close()
H.close()

#E=EllipticCurve()
