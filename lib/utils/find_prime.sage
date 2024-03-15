from argparse import ArgumentParser
parser=ArgumentParser()

parser.add_argument('-c','--conductor')
parser.add_argument('-n','--nbprimes')

args=parser.parse_args()

conductor=Integer(args.conductor)
nb_primes=Integer(args.nbprimes)


filename="../candidate_conductors/candidate_conductors"+str(nb_primes)+".md"
g=open(filename, "r")
line=0

nim=abs(Integer(g.readline().split(" ")[0]))
nre=abs(Integer(g.readline().split(" ")[0]))

while nim!=conductor and nre!=conductor:
	nim=abs(Integer(g.readline().split(" ")[0]))
	nre=abs(Integer(g.readline().split(" ")[0]))
	line+=2

print(line)
print(conductor,nim, nre)
