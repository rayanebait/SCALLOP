K.<t>=GF(49)
E0=EllipticCurve(K, [-1,0])

print(E0)

P=E0.lift_x(1+t)
Q=E0.lift_x(4)
print(P.order(), Q.order())
NPS=[n*P for n in range(1, 15)]
print(NPS)

