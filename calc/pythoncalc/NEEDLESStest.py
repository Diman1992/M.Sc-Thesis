import numpy as np
from sympy import *

l = 0.225
A = 0.814
r = 0.135
n = 0.349j

print(r-n)
ckm = np.zeros((3,3), dtype=np.complex_)
v = np.zeros((3,3))
yuk = np.zeros((3,3))

ckm[0][0] = 1-l**2/2
ckm[0][1] = l
ckm[0][2] = A*l**3*(r-n)
ckm[1][0] = -l
ckm[1][1] = 1-l**2/2
ckm[1][2] = A*l**2
ckm[2][0] = A*l**3*(1-r-n)
ckm[2][1] = -A*l**2
ckm[2][2] = 1

v[0][0] = 1
v[0][1] = l
v[0][2] = l**3
v[1][0] = -l
v[1][1] = 1
v[1][2] = l**2
v[2][0] = -l**3
v[2][1] = -l**2
v[2][2] = 1

e = 0.2
qi = np.array([3,2,0])
ui = np.array([3,2,0])
for i in range(0,3):
	for j in range(0,3):
		yuk[i][j] = e**(qi[i]+ui[j])
yukTrans = np.transpose(yuk)
res = np.matrix.dot(yuk,yukTrans)

print(yuk)
print(yukTrans)
print(res)
det = np.linalg.det(yuk)
print(det)
