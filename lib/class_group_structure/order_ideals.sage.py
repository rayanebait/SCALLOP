

# This file was *autogenerated* from the file order_ideals.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_96754657397296242025209609489487721053638409 = Integer(96754657397296242025209609489487721053638409)
K = QuadraticField(-_sage_const_1 , names=('t',)); (t,) = K._first_ngens(1)

#62516645812060079426901900471388824796582670407138252131237014712949
O=K.order(_sage_const_96754657397296242025209609489487721053638409 *t)

print(O.class_number())


