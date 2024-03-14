def FullRepresentInteger(M, p, trials=0):
	m_=floor(float(sqrt(4*M/p)))
	z_=random.randint(-m_, m_)
	m__=floor(float( sqrt(4*M/p-z_**2) ))
	t_=random.randint(-m__, m__)
	M_=4*M-p*((z_)**2+(t_)**2)
	try:
		x_,y_=two_squares(M_)
	except ValueError:
		if trials>1000:
			return
		FullRepresentInteger(M,p,trials+1)
	if (x_-t_)%2!=0 or (z_-y_)%2!=0:
		if trials>1000:
			return
		FullRepresentInteger(M,p,trials+1)
	return (x_, y_, z_, t_)


