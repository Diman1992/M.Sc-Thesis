import numpy as np
from sympy import *

u = np.zeros((3,3))
s12sq = 0.312

s12 = np.sqrt(s12sq)
c12 = np.sqrt(1 - s12sq)
sqr2 = np.sqrt(2)

u[0][0] = c12
u[0][1] = s12
u[0][2] = 0
u[1][0] = -s12/sqr2
u[1][1] = c12/sqr2
u[1][2] = -1/sqr2
u[2][0] = -s12/sqr2
u[2][1] = c12/sqr2
u[2][2] = 1/sqr2
udag = np.transpose(u)
#sympy
a1 = Symbol('x')
a2 = Symbol('y')
a3 = Symbol('y')
a4 = Symbol('y')
a5 = Symbol('z')
a6 = Symbol('w')
a7 = Symbol('y')
a8 = Symbol('w')
a9 = Symbol('z')

u1 = Symbol('c12')
u2 = Symbol('s12')
u3 = Symbol('0')
u4 = Symbol('-s12/q2')
u5 = Symbol('c12/q2')
u6 = Symbol('-1/q2')
u7 = Symbol('-s12/q2')
u8 = Symbol('c12/q2')
u9 = Symbol('1/q2')

##u2 = s12
#u3 = 0
#u4 = -s12/sqr2
#u5 = -c12/sqr2
#u6 = -1/sqr2
#u7 = -s12/sqr2
#u8 = c12/sqr2
#u9 = 1/sqr2


#a = Matrix([[a1,a2,a3],[a4,a5,a6],[a7,a8,a9]])
#res = a*u
#res = udag.dot(res)
#u1 = c12
#u = Matrix([[u1,u2,u3],[u4,u5,u6],[u7,u8,u9]])
#udag = Matrix([[u1,u4,u7],[u2,u5,u8],[u3,u6,u9]])
#res =udag* a*u
#print(res)


######S_TB
s = np.zeros((3,3))
a = np.zeros((3,3))
for i in range(0,3):
	for j in range(0,3):
		s[i][j] = 2
		if(i==j):
			s[i][j] = -1
		a[i][j] = 2*i+j+1
s = 0.333*s

b = s.dot(a.dot(s))

print(a)	
print(b)
