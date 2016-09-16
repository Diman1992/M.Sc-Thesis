import numpy as np
import random

def makeA(randOr1): #coefficients of the entries of the yukawa-matrix
	np.random.seed()
	aij = np.zeros((3,3),dtype=complex)
	if(randOr1==0):
		for i in range(0,3):
			for j in range(0,3):
				r, phi = np.random.uniform(0.9,1.1), np.random.uniform(0,2*np.pi)
				aij[i][j] = r*np.exp(phi*1j)
	if(randOr1==1):
		for i in range(0,3):
			for j in range(0,3):
				aij[i][j] = np.random.uniform(0.9,1.1)
	if(randOr1==2):
		for i in range(0,3):
			for j in range(0,3):
				aij[i][j] = 1
	return aij
	
def makeN(di, si): #exponents of the yukawa entries
	nExp = np.zeros((3,3))
	for i in range(0,3):
		for j in range(0,3):
			nExp[i][j] = qi[i] + si[j]	
	return nExp

def makeYs(aij, di, si, eps): #producing yukawa matricies by coefficients a and charges n
	nExp = makeN(qi,ui)
	Y = aij*eps**nExp 
	return Y

def makeVs_NumUnitary(Y): #the rotation matricies should be unitary, here computed by the recipe of the paper
	V_L, V_R = np.zeros((3,3),dtype=complex), np.zeros((3,3),dtype=complex)
	for i in range(0,3):
		for j in range(0,3):
			if(j>=i):
				V_L[i][j] = Y[i][j]/Y[j][j]
				V_R[i][j] = Y[j][i]/Y[j][j]
	#unitarity conditions
	V_L[1][0] = -np.conjugate(V_L[0][1]) - V_L[1][2]*np.conjugate(V_L[0][2])
	V_L[2][0] = (np.conjugate(V_L[1][2])*np.conjugate(V_L[0][1])+np.conjugate(V_L[0][2])) / (1-np.conjugate(V_L[1][0]*V_L[0][1]))
	V_L[2][1] = -V_L[2][0]*np.conjugate(V_L[1][0]) - np.conjugate(V_L[1][2])

	V_R[1][0] = -np.conjugate(V_R[0][1]) - V_R[1][2]*np.conjugate(V_R[0][2])
	V_R[2][0] = (np.conjugate(V_R[1][2])*np.conjugate(V_R[0][1])+np.conjugate(V_R[0][2])) / (1-np.conjugate(V_R[1][0]*V_R[0][1]))
	V_R[2][1] = -V_R[2][0]*np.conjugate(V_R[1][0]) - np.conjugate(V_R[1][2])

#	V_L, V_R = normaliseVCols(V_L), normaliseVCols(V_R)
	return V_L, V_R

def makeVs_AnaBidiag(Y): #the rotation matricies should be unitary, here computed analytically by bidiagonalisation (U* MM* U = M2diag)
	V_L, V_R = np.zeros((3,3),dtype=complex), np.zeros((3,3),dtype=complex)
	Ynt = Y.dot(Y.transpose().conjugate()) #used for V_L
	Ytn = Y.transpose().conjugate().dot(Y) #used for V_R
	evYnt, evYtn = np.linalg.eig(Ynt)[1], np.linalg.eig(Ytn)[1]
	for i in range(0,3): #swap columns to get the O(1) entries at the diagonal
		V_L[i][0] = evYnt[i][1]
		V_L[i][1] = evYnt[i][2]
		V_L[i][2] = evYnt[i][0]
		V_R[i][0] = evYtn[i][1]
		V_R[i][1] = evYtn[i][2]
		V_R[i][2] = evYtn[i][0]
	return V_L, V_R

def make_AnaDiag(m):
	P = np.zeros((3,3))
	ev = np.linalg.eig(m)[1]
	for i in range(0,3):
		P[i][0] = ev[i][1]
		P[i][1] = ev[i][2]
		P[i][2] = ev[i][0]
	return P

def checkDiag(V_L, Y, V_R):
	return np.transpose(V_L).dot(Y).dot(V_R) #should return an (almost) diagonal yukawa matrix

def det(V): #determinant to prove another property of unitary matricies
	print(np.absolute(np.linalg.det(V)))

def normaliseVCols(V): #in case of degenerate charges (as for Vd_R here) the columns have to be normalised so that V becomes unitary
	norm = np.zeros((3))
	for j in range(0,3):
		for i in range(0,3):
			norm[j] += V[i][j]**2
	norm = np.sqrt(norm)
	for j in range(0,3):	
		for i in range(0,3):
			V[i][j] = 1/norm[j] * V[i][j]
	print(norm)
	return V
	
def printVs(Vu_L, Vu_R, Vd_L, Vd_R, Vl_L, Vl_R, nOra):
	print("#######Rotmatrices######")
	print("######%s#######") % (nOra)
	print(np.absolute(Vu_L))
	print(np.absolute(Vu_R))
	print(np.absolute(Vd_L))
	print(np.absolute(Vd_R))
	print(np.absolute(Vl_L))
	print(np.absolute(Vl_R))
	print("************************")

def printYnorm(Yu, Yd, Yl, nOra):
	print("#######Ynorm######")
	print("######%s#######") % (nOra)
	print(np.absolute(Yu))
	print(np.absolute(Yd))
	print(np.absolute(Yl))
	print("************************")

def printVdet(vul,vur,vdl,vdr,vll,vlr,nOra):
	print("######Determinants####")
	print("######%s#######") % (nOra)
	det(vul),det(vur),det(vdl),det(vdr),det(vll),det(vlr)
	print("************************")

def func(A_ass, eps, qi, ui, di, li, ei): #main function
	aiju, aijd, aijl = makeA(A_ass), makeA(A_ass), makeA(A_ass)
	Yu, Yd, Yl = makeYs(aiju,qi,ui,eps),makeYs(aijd,qi,di,eps),makeYs(aijl,li,ei,eps)

	####paper####
	Vu_Lnum, Vu_Rnum = makeVs_NumUnitary(Yu)
	Vd_Lnum, Vd_Rnum = makeVs_NumUnitary(Yd)
	Vl_Lnum, Vl_Rnum = makeVs_NumUnitary(Yl)
	YuNorm_num = checkDiag(Vu_Lnum, Yu, Vu_Rnum)
	YdNorm_num = checkDiag(Vd_Lnum, Yd, Vd_Rnum)
	YlNorm_num = checkDiag(Vl_Lnum, Yl, Vl_Rnum)
		
	####mathematics####
	Vu_Lana, Vu_Rana = makeVs_AnaBidiag(Yu)
	Vd_Lana, Vd_Rana = makeVs_AnaBidiag(Yd)
	Vl_Lana, Vl_Rana = makeVs_AnaBidiag(Yl)
	YuNorm_ana = checkDiag(Vu_Lana, Yu, Vu_Rana)
	YdNorm_ana = checkDiag(Vd_Lana, Yd, Vd_Rana)
	YlNorm_ana = checkDiag(Vl_Lana, Yl, Vl_Rana)

	####printresults####
	print("always uL, uR, dL, dR")
	printVdet(Vu_Lana,Vu_Rana,Vd_Lana,Vd_Rana,Vl_Lana,Vl_Lana,"ana")
	printVdet(Vu_Lnum,Vu_Rnum,Vd_Lnum,Vd_Rnum,Vl_Lnum,Vl_Rnum,'num')
	printVs(Vu_Lana,Vu_Rana,Vd_Lana,Vd_Rana,Vl_Lana,Vl_Rana,'ana')
	printVs(Vu_Lnum,Vu_Rnum,Vd_Lnum,Vd_Rnum,Vl_Lnum,Vl_Rnum,'num')
	printYnorm(YuNorm_ana, YdNorm_ana,YlNorm_ana, 'ana')
	printYnorm(YuNorm_num, YdNorm_num,YlNorm_num, 'num')

##############################
######## assignments #########
##############################
qi = np.array([3,2,0]) #charges under the U(1)_F symmetry for quark doublets
ui = np.array([3,2,0]) #up-type singlets
di = np.array([4,2,2]) #down-type singlets
li = np.array([3,1,0]) #lepton doublets
ei = np.array([5,4,3]) #charged lepton singlets

YukPrefactor = 0 #0complex(generic), 1real, 2unity
expansion = 0.2 #Cabibbo parameter

######
xi, m = np.zeros((3,3)), np.zeros((3,3))
m[0][0] = 8
m[0][1] = 5
m[0][2] = 6
m[1][0] = 5
m[1][1] = 4
m[1][2] = 4
m[2][0] = 6
m[2][1] = 4
m[2][2] = 2
m = expansion**m

xi[0][1] =1
xi[1][0] =1


print(m)
p = make_AnaDiag(m)
Vl, Vr = makeVs_AnaBidiag(m)
print(np.linalg.eig(m)[0])
print(Vl)
print(Vr)
print(checkDiag(Vl, m, Vr))
print(xi)
print(Vl.dot(xi).dot(Vr.transpose()))
print(Vl.transpose().dot(xi).dot(Vr))





#func(YukPrefactor, expansion, qi, ui, di, li, ei)  

#EOP
