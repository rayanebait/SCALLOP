from multiprocessing import Process, TimeoutError, Pipe

import time

def f(cqr, fcs):
	while True:
		print("Waiting")
		i=cqr.recv()

		print("Factoring")
		factors=ecm.factor(11111111111111111111111111111111111111111111111111111111111111111111+2*i, None)

		fqs.send(factors)


processes=[]


cqr, cqs=Pipe(duplex=False)
fqr, fqs=Pipe(duplex=False)

p=Process(target=f, args=(cqr,fqs))

p.start()

i=0
while True:
	time.sleep(float(0.100))
	cqs.send(i)

	i+=1
	try:
		time.sleep(float(0.2))
		if fqr.poll():
			factors=fqr.recv()
		else:
			p.terminate()
			print("process terminated\n")
			cqr, cqs=Pipe(duplex=False)
			fqr, fqs=Pipe(duplex=False)

			p=Process(target=f, args=(cqr,fqs))
			p.start()
			continue

	except TimeoutError:
		p.terminate()
		print("process terminated\n")
		cqr, cqs=Pipe(duplex=False)
		fqr, fqs=Pipe(duplex=False)

		p=Process(target=f, args=(cqr,fqs))
		p.start()
		continue
	print("has factored",factors)

	



#time.sleep(1)
print("DONE")
