from fpylll import * 
from copy import copy
n=10
m=10

FPLLL.set_random_seed(10932)
#IntegerMatrix.random.__code__.co_varnames
A=IntegerMatrix.random(10, "a", k=10, bits=100)

B = BKZ.reduction(copy(A), BKZ.EasyParam(40, max_loops=8))
print(A)
