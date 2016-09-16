import numpy as np

def printMat(mat, name):
	print("####%s####") % (name)
	print(mat)

def howUnitaryIsMat(mat,dim):
	norms = list()
	for j in range(0,dim):
		for k in range(j,dim):
			l = list()
			n = 0
			for i in range(0,dim):
				n += mat[i][j]*mat[i][k]
			l.append(j)
			l.append(k)
			l.append(n)
			norms.append(l)
	one = mat.dot(mat.transpose())
	print(mat)
	print(norms)
	print(one)
	print(np.linalg.det(mat))

i = 0
if (i ==0): #S3
	s1 = np.zeros((3,3),dtype=complex)
	s2 = np.zeros((3,3),dtype=complex)
	s3 = np.zeros((3,3),dtype=complex)
	
	s1[0][0] = 0
	s1[0][1] = 1
	s1[0][2] = 0
	s1[1][0] = 1
	s1[1][1] = 0
	s1[1][2] = 1
	s1[2][0] = 0
	s1[2][1] = 1
	s1[2][2] = 0
	
	s2[0][0] = 0
	s2[0][1] = -1j
	s2[0][2] = 0
	s2[1][0] = 1j
	s2[1][1] = 0
	s2[1][2] = -1j
	s2[2][0] = 0
	s2[2][1] = 1j
	s2[2][2] = 0

	s3[0][0] = 1
	s3[0][1] = 0
	s3[0][2] = 0
	s3[1][0] = 0
	s3[1][1] = 0
	s3[1][2] = 0
	s3[2][0] = 0
	s3[2][1] = 0
	s3[2][2] = -1

	print(s1.dot(s2)-s2.dot(s1))
	print(s2.dot(s3)-s3.dot(s2))
	print(s3.dot(s1)-s1.dot(s3))

if (i==1): #T'
	s = np.zeros((3,3),dtype=complex)
	t = np.zeros((3,3),dtype=complex)
	w = np.exp(2*np.pi/3 *1j)
	s[0][0]	= -1
	s[0][1]	= 2*w
	s[0][2]	= 2*w**2
	s[1][0]	= 2*w**2
	s[1][1]	= -1
	s[1][2]	= 2*w
	s[2][0]	= 2*w
	s[2][1]	= 2*w**2
	s[2][2]	= -1
	t[0][0]	= 1
	t[0][1]	= 0
	t[0][2]	= 0
	t[1][0]	= 0
	t[1][1]	= w
	t[1][2]	= 0
	t[2][0]	= 0
	t[2][1]	= 0
	t[2][2]	= w**2
	s = 0.333*s

	m = np.zeros((3,3),dtype=complex)
	m[1][0] = 1
	m[1][1] = 1
	m[1][2] = 1
	m[2][1] = 1
	m[2][2] = 1

	print(t.dot(m).dot(t.transpose()))
if (i==2): #majorana spinor
	g0, g1, g2, g3, g5, u = np.zeros((4,4),dtype=complex),np.zeros((4,4),dtype=complex),np.zeros((4,4),dtype=complex),np.zeros((4,4),dtype=complex),np.zeros((4,4),dtype=complex),np.zeros((4,4),dtype=complex),
	g0[0][2] = 1
	g0[1][3] = 1
	g0[2][0] = 1
	g0[3][1] = 1
	g1[0][3] = 1
	g1[1][2] = 1
	g1[2][1] = -1
	g1[3][0] = -1
	g2[0][3] = -1j
	g2[1][2] = 1j
	g2[2][1] = 1j
	g2[3][0] = -1j
	g3[0][2] = 1
	g3[1][3] = -1
	g3[2][0] = -1
	g3[3][1] = 1
	g5[0][0] = -1
	g5[1][1] = -1
	g5[2][2] = 1
	g5[3][3] = 1
	u[0][2] = 1
	u[1][3] = 1
	u[2][0] = -1
	u[3][1] = -1
	m = list()
	m.append(g0)
	m.append(g1)
	m.append(g2)
	m.append(g3)
	m.append(g5)
#	m.append(u)
	a = np.array((4,4),dtype=complex)
	for i in range(0, len(m)):
		for j in range(0,len(m)):
			for k in range(0,len(m)):
				for l in range (0,len(m)):
					for q in range(0,len(m)):
						a = m[i].dot(m[j]).dot(m[k]).dot(m[l]).dot(m[q])	
						if((1j*a == u).all()):	
								print(i,j,k,l,q)	
	print(1j*(g1.dot(g2).dot(g3)))

if(i==3):	#massdiag
	e1 = Symbol('(y1T2)')
	e2 = Symbol('(y1T3+y2u)')
	e3 = Symbol('(-y4v2)')
	e4 = Symbol('(y1T3-y2u)')
	e5 = Symbol('(y1T1)')
	e6 = Symbol('(y4v1)')
	e7 = Symbol('(-y3v2)')
	e8 = Symbol('(y3v1)')
	e9 = Symbol('(ytL)')
	e1 = 0
	e9 = 1

	m = Matrix(3,3,[e1,e2,e3,e4,e5,e6,e7,e8,e9])
	x = m.eigenvals()
	print("x",x)

if(i==4):
	m = np.ones((3,3))
	l = 0.2
	m[0][0] = -l
	m[0][2] = -l**2
	m[1][1] = l
	m[1][2] = -l**3
	m[2][0] = -l**2
	m[2][1] = -l**3
	m[2][2] = l**5
	howUnitaryIsMat(m,3)

if(i==5):
	u = np.zeros((3,3))
	u[0][0] = np.sqrt(2.0/3)
	u[0][1] = 1/np.sqrt(3)
	u[0][2] = 0
	u[1][0] = -1/np.sqrt(6)
	u[1][1] = 1/np.sqrt(3)
	u[1][2] = -1/np.sqrt(2)
	u[2][0] = -1/np.sqrt(6)
	u[2][1] = 1/np.sqrt(3)
	u[2][2] = 1/np.sqrt(2)

	g = u.dot(u.transpose())
	print(g)
