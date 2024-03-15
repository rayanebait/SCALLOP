from argparse import ArgumentParser
parser=ArgumentParser()

parser.add_argument('-v','--verbose', action='store_true')
parser.add_argument('-c','--conductor')
parser.add_argument('-n','--nbprimes')

args=parser.parse_args()

conductor=Integer(args.conductor)
nb_primes=Integer(args.nbprimes)


filename="../candidate_conductors/candidate_conductors"+str(nb_primes)+".md"
g=open(filename, "r")
line=0

if args.verbose:
	print(f"Searching conductor {conductor} in file {filename}\n")

nim=abs(Integer(g.readline().split(" ")[0]))
nre=abs(Integer(g.readline().split(" ")[0]))

while nim!=conductor and nre!=conductor:
	nim=abs(Integer(g.readline().split(" ")[0]))
	nre=abs(Integer(g.readline().split(" ")[0]))
	line+=2

if args.verbose:
	print(f"Found conductor at line {line} with:\n\t {nre}\n\t {nim}\n")
else:
	print(conductor,"\n",nim,"\n", nre,"\n")
