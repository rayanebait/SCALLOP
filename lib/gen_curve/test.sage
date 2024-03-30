from argparse import ArgumentParser

parser=ArgumentParser()

parser.add_argument("-n", "--regularity", default="3")
parser.add_argument("-d", "--depth", default="10")
parser.add_argument("-i", "--iterations", default="100")



args=parser.parse_args()

def graphs_n_reg_depth_h(n,h, base_n_int, step=0, iterations=0):
	def do_stuff():
		print(base_n_int)
	
	i=0
	while i<iterations:
		print(step)
		while step < h-1:
			step+=1
	
		do_stuff()
		while base_n_int[step]<n-1:
			base_n_int[step]+=1
			do_stuff()

		while base_n_int[step]==n-1:
			base_n_int[step]=0
			step-=1
			if base_n_int[step]==n-1:
				continue
			base_n_int[step]+=1
		i+=3


graphs_n_reg_depth_h(Integer(args.regularity),\
			Integer(args.depth),\
			[0]*Integer(args.depth),\
			iterations=Integer(args.iterations))

