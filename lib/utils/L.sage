def L(x,a=1/2,c=1): return exp(c*( (ln(x)**a) * (ln(ln(x))**(1-a)) ))

log_2=80

value=2*floor(float(L(2**log_2)))

f=open("../Lbound.md", 'w', encoding="utf-8")
f.write(f'{value}\n')

f.close()
