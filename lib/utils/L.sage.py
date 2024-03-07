

# This file was *autogenerated* from the file ../utils/L.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_80 = Integer(80)
def L(x,a=_sage_const_1 /_sage_const_2 ,c=_sage_const_1 ): return exp(c*( (ln(x)**a) * (ln(ln(x))**(_sage_const_1 -a)) ))

log_2=_sage_const_80 

value=_sage_const_2 *floor(float(L(_sage_const_2 **log_2)))

f=open("../Lbound.md", 'w', encoding="utf-8")
f.write(f'{value}\n')

f.close()
