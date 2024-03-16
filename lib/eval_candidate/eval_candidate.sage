from multiprocessing import Process, Pipe 
import time

def LBound(x,a=1/2,c=1): return exp(c*( (ln(x)**a) * (ln(ln(x))**(1-a)) ))

def is_L_smooth(factors, L):
	is_smooth=True
	for i in factors:
		if i>10*L:
			is_smooth=False
			break
	return is_smooth


def fac(cpr,fps):
	while True:
		#print("Waiting")
		c=cpr.recv()
		c=abs(c)
		L=floor(float(LBound(c)))
		#print("Factoring")

		if c&3==1:
			factors=ecm.factor(c-1,None)
		else:
			factors=ecm.factor(c+1,None)

		fps.send((factors, c, L))


def eval_candidates_for_nprimes(nprimes):
	name="../candidate_conductors/candidate_conductors"+str(nprimes)+".md"
	F=open(name, "r", encoding="utf-8")

	nprimes=Integer(F.readline())
	nb_factored=0
	nb_candidates=15

	candidates=[]
	fs=[2]*(2*nb_candidates)

	f_re_str=''
	f_im_str=''

	processes=[None]*(2*nb_candidates)

	pipes_cand=[]
	pipes_fac=[]

	#cqr, cqs

	for i in range(2*nb_candidates):
		pipes_cand.append(Pipe(duplex=False))
		pipes_fac.append(Pipe(duplex=False))
		processes[i]=Process(target=fac, args=(pipes_cand[i][0],pipes_fac[i][1]))
		processes[i].start()


	line=0
	eof=False
	while not eof:
		line+=2*nb_candidates
		nb_factored+=1
		#Multiple loops for readability, doesn't impact perf here

		for i in range(nb_candidates):
		#Read nb_processes candidates
			f_re_str=F.readline()
			f_im_str=F.readline()
			if f_re_str == '' or f_im_str=='':
				eof=True
				break
			f_re_str=f_re_str.split(' ')
			f_im_str=f_im_str.split(' ')
			
			fs[i]=(Integer(f_re_str[0]),Integer(f_re_str[1][0]))
			fs[i+1]=(Integer(f_im_str[0]),Integer(f_im_str[1][0]))
			print("Attempting to factor:",fs[i], fs[i+1])
		if eof:
			break

		time.sleep(float(0.1))
		for i in range(nb_candidates):
		#send each candidate to a process
			if fs[i][1]==0:
				continue
			pipes_cand[i][1].send(fs[i][0])

		time.sleep(float(0.3))
		for i in range(nb_candidates):
		#Check only if candidate was prime
			if fs[i][1]==0:
				continue
			#timeout at 0.5 secs
			if pipes_fac[i][0].poll():
				#check if process nb i terminated
				factors=pipes_fac[i][0].recv()
			else:
				#terminate the process and launch another one
				processes[i].terminate()
				print(f"Candidate {fs[i]} timed out\n")

				pipes_cand[i]=Pipe(duplex=False)
				pipes_fac[i]=Pipe(duplex=False)
				processes[i]=Process(target=fac, args=(pipes_cand[i][0],pipes_fac[i][1]))
				processes[i].start()

				continue

			#check if candidate has L-smooth class number
			if is_L_smooth(factors[0], factors[2]):
				print(f"Conductor {factors[1]} gives: {factors[0]}")
				candidates.append(factors)
				if len(candidates)==5:
					eof=True
					break
			else:
				print(f"Candidate {factors[1]} doesnt give {factors[2]}-smooth class group:\n\t{factors[0]}")

	time.sleep(float(0.100))
	for i in range(len(fs)):
	#send each candidate to a process
		if type(fs[i])==Integer or fs[i][1]==0:
			continue
		pipes_cand[i][1].send(fs[i][0])

	time.sleep(float(0.8))
	for i in range(len(fs)):
		#timeout at 0.3 secs
		if type(fs[i])==Integer or fs[i][1]==0:
			continue
		if pipes_fac[i][0].poll():
			#check if process nb i terminated
			factors=pipes_fac[i][0].recv()
		else:
			#terminate the process and launch another one
			processes[i].terminate()
			print(f"Candidate {fs[i]} timed out\n")

			pipes_cand[i]=Pipe(duplex=False)
			pipes_fac[i]=Pipe(duplex=False)
			processes[i]=Process(target=fac, args=(pipes_cand[i][0],pipes_fac[i][1]))
			processes[i].start()
			continue

		#check if candidate has L-smooth class number
		if is_L_smooth(factors[0], factors[2]):
			print(f"Conductor {factors[1]} gives: {factors[0]}")
			candidates.append(factors)
			#if len(candidates)==3:
				#eof=True
				#break
		else:
			print(f"Candidate {factors[1]} doesnt give {factors[2]}-smooth class group:\n\t{factors[0]}")


	for p in processes:
		p.terminate()
	F.close()
	print("Number of factored candidates", nb_factored*nb_candidates)
	print(candidates)

	out_file_name="conductors_data"+str(nprimes)+".md"
	G=open(out_file_name, "w", encoding="utf-8")
	G.write(f'Number of primes used: {nprimes}\n')
	for candidate in candidates:
		flog=floor(float(log(candidate[1])/log(2)))
		faclog=floor(float(log(candidate[0][-1])/log(2)))
		Llog=floor(float(log(candidate[2])/log(2)))

		G.write(f"{candidate[1]} is O({candidate[2]})-smooth:\n\t {candidate[0]}\n")
		G.write(f"{candidate[1]} has size 2**{flog} and 2**{faclog}-smooth with L(f,1/2) about 2**{Llog}\n\n")
	G.write(f"finished at line: {line}\n")
	
	G.close()

nprimes=26

while nprimes < 27:
	eval_candidates_for_nprimes(nprimes)
	nprimes+=1
