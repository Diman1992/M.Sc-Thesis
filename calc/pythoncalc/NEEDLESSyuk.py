import numpy as np
from sympy import *
import random
e = 0.2
### with numpy

def sympy(): #function=0
	
	e1 = Symbol('e11')
	e2 = Symbol('e12')
	e3 = Symbol('e13')
	e4 = Symbol('e21')
	e5 = Symbol('e22')
	e6 = Symbol('e23')
	e7 = Symbol('e31')
	e8 = Symbol('e32')
	e9 = Symbol('e33')

	e1 = 1
#	e2 = e**(1)
	e5 = 1
#	e3 = e**(3)
#	e6 = e**(2)
	e9 = 1

	test = Matrix([[e1,e2,e3],[e4,e5,e6],[e7,e8,e9]])
	testHerm = test.H
	print("#########################################")
	print("############sympy########################")
	print("#########################################")
	print("Matrix: ")
	print(test)
	print("HermMat: ")
	print(testHerm)
	print("###############################")
	print("Result: ")
	print(test.dot(testHerm))
	print("ResultSwapped: " ) #for commutationproperties
	print(testHerm.dot(test))
	print("Determinant:")
	print(test.det())
	print("*****************************************")
	return test, testHerm

def diagonality (upOrDown): #function=1, down=0, up=1
	### check diagonal Y = V_L Y V_R
	#Flavourcharges
	qi = np.array([3,2,0]) 
	ui = np.array([3,2,0])
	di = np.array([4,2,2])

	nExpU = np.zeros((3,3))
	nExpD = np.zeros((3,3))
	#nExpU are the indices in the eps-exponent for Y,V_L,V_R
	for i in range(0,3):
		for j in range(0,3):
			nExpU[i][j] = qi[i] + ui[j]
			nExpD[i][j] = qi[i] + di[j]

	aiju, aijd = createA(1)[0], createA(1)[1]
	#Yukawa matrix with flavourcharges
	Yu = aiju*np.array(e**nExpU)
	Yd = aijd*np.array(e**nExpD)

	#Rotationmatrices V_L and V_R for up- and downtype
	Vu_L = np.zeros((3,3))
	Vd_L = np.zeros((3,3))
	Vu_R = np.zeros((3,3))
	Vd_R = np.zeros((3,3))
	for i in range(0,3):
		for j in range(0,3):
			if(j>=i):
				Vu_L[i][j] = Yu[i][j]/Yu[j][j]
				Vd_L[i][j] = Yd[i][j]/Yd[j][j]
				Vu_R[i][j] = Yu[j][i]/Yu[j][j]
				Vd_R[i][j] = Yd[j][i]/Yd[j][j]
			
	#making Vs unitary (transpose because j<->i and hardcoded (following) unitarity formula was already made, so theyre not getting overwritten)
	print("first triangle:",Vu_L)
	Vu_L=Vu_L.transpose()
	Vd_L=Vd_L.transpose()
	Vu_R=Vu_R.transpose()
	Vd_R=Vd_R.transpose()

	Vu_L[0][1] = -(Vu_L[1][0]*Vu_L[1][1]+Vu_L[2][0]*Vu_L[2][1])/(Vu_L[0][0])
	Vu_L[1][2] = ((Vu_L[0][1]*Vu_L[2][1]*Vu_L[2][2])/(Vu_L[0][0]) - (Vu_L[2][1]*Vu_L[2][2])/(Vu_L[1][1])) / ((1-(Vu_L[0][1]*Vu_L[1][0])/(Vu_L[1][1]*Vu_L[0][0])))
	Vu_L[0][2] = -(Vu_L[1][0]*Vu_L[1][2]+Vu_L[2][0]*Vu_L[2][2])/(Vu_L[0][0])

	Vd_L[0][1] = -(Vu_L[1][0]*Vu_L[1][1]+Vu_L[2][0]*Vu_L[2][1])/(Vu_L[0][0])
	Vd_L[1][2] = ((Vu_L[0][1]*Vu_L[2][1]*Vu_L[2][2])/(Vu_L[0][0]) - (Vu_L[2][1]*Vu_L[2][2])/(Vu_L[1][1])) / ((1-(Vu_L[0][1]*Vu_L[1][0])/(Vu_L[1][1]*Vu_L[0][0])))
	Vd_L[0][2] = -(Vu_L[1][0]*Vu_L[1][2]+Vu_L[2][0]*Vu_L[2][2])/(Vu_L[0][0])


	Vu_R[0][1] = -(Vu_R[1][0]*Vu_R[1][1]+Vu_R[2][0]*Vu_R[2][1])/(Vu_R[0][0])
	Vu_R[1][2] = ((Vu_R[0][1]*Vu_R[2][1]*Vu_R[2][2])/(Vu_R[0][0]) - (Vu_R[2][1]*Vu_R[2][2])/(Vu_R[1][1])) / ((1-(Vu_R[0][1]*Vu_R[1][0])/(Vu_R[1][1]*Vu_R[0][0])))
	Vu_R[0][2] = -(Vu_R[1][0]*Vu_R[1][2]+Vu_R[2][0]*Vu_R[2][2])/(Vu_R[0][0])

	Vd_R[0][1] = -(Vd_R[1][0]*Vd_R[1][1]+Vd_R[2][0]*Vd_R[2][1])/(Vd_R[0][0])
	Vd_R[1][2] = ((Vd_R[0][1]*Vd_R[2][1]*Vd_R[2][2])/(Vd_R[0][0]) - (Vd_R[2][1]*Vd_R[2][2])/(Vd_R[1][1])) / ((1-(Vd_R[0][1]*Vd_R[1][0])/(Vd_R[1][1]*Vd_R[0][0])))
	Vd_R[0][2] = -(Vd_R[1][0]*Vd_R[1][2]+Vd_R[2][0]*Vd_R[2][2])/(Vd_R[0][0])

#	Vd_R = 0.5**(0.333) *  Vd_R #getting abs(det) = 1
	print("after hardcoding unitarity:",Vu_L)
	#transposition 
	Vu_Ltrans = Vu_L.transpose()
	Vd_Ltrans = Vd_L.transpose()

	#matrixmultiplikation
	YuNorm = np.dot(np.dot(Vu_Ltrans, Yu),Vu_R)
	YdNorm = np.dot(np.dot(Vd_Ltrans, Yd),Vd_R)
	if(upOrDown == 0):
		print"##########################################"
		print"#######Calculation for Yukawa#############"
		print"V_L+"
		print(Vd_Ltrans)
		print("Yd")
		print(Yd)
		print("V_R")
		print(Vd_R)
		print("eigenvalues of Yd", np.linalg.eig(Yd)[0])
		print"###############Result#####################"
		print"Ynorm"
		print(YdNorm)
		return Vd_L, Vd_Ltrans, Vd_R, Vd_R.transpose()
	if(upOrDown == 1):
		print"##########################################"
		print"#######Calculation for Yukawa#############"
		print"V_L+"
		print(Vu_Ltrans)
		print("Yu")
		print(Yu)
		print("V_R")
		print(Vu_R)
		print("eigenvalues of Yu", np.linalg.eig(Yu)[0])
		print"###############Result#####################"
		print"Ynorm"
		print(YuNorm)
		return Vu_L, Vu_Ltrans, Vu_R, Vu_R.transpose()
def unitarityOfV(param,lOrR):
	V = diagonality(param)[lOrR] 	
	Vdagger = diagonality(param)[lOrR+1]
	print("#####unitarityCheck", lOrR,"#####")
	print("V: ")
	print(V)
	print("Vdag: ")
	print(Vdagger)
	print("UnitRes: ")
	print(Vdagger.dot(V))
	print("Determinant: ")
	print(np.linalg.det(V))

def createA(oneOrRandom):
	np.random.seed(1)
	aiju, aijd = np.zeros((3,3)), np.zeros((3,3))
	if(oneOrRandom==0):
		for i in range(0,3):
			for j in range(0,3):
				aiju[i][j] = np.random.uniform(0.8,1.2)
				aijd[i][j] = np.random.uniform(0.8,1.2)
	if(oneOrRandom==1):
		for i in range(0,3):
			for j in range(0,3):
				aiju[i][j] = 1
				aijd[i][j] = 1
	return aiju, aijd
		
	
def biUnitaryDiag(): #function=4
	qi = np.array([3,2,0])
	ui = np.array([3,2,0])
	di = np.array([4,2,2])
	
	nExpU = np.zeros((3,3))
	nExpD = np.zeros((3,3))
	aiju, aijd = createA(1)[0], createA(1)[0]
	for i in range(0,3):
		for j in range(0,3):
			nExpU[i][j] = qi[i] + ui[j]
			nExpD[i][j] = qi[i] + di[j]
		
	Yu, Yd = aiju*np.array(e**nExpU), aijd*np.array(e**nExpD)


	Yu2normaltrans = Yu.dot(Yu.transpose()) #for VuL
	Yu2transnormal = Yu.transpose().dot(Yu) #for VuR
	Yd2normaltrans = Yd.dot(Yd.transpose()) #for VdL
	Yd2transnormal = Yd.transpose().dot(Yd) #for VdR
	
	evYu2nt, evYu2tn = np.linalg.eig(Yu2normaltrans)[1], np.linalg.eig(Yu2transnormal)[1]
	evYd2nt, evYd2tn = np.linalg.eig(Yd2normaltrans)[1], np.linalg.eig(Yd2transnormal)[1]

	VuL = np.zeros((3,3))
	VuR = np.zeros((3,3))
	VdL = np.zeros((3,3))
	VdR = np.zeros((3,3))
	for i in range(0,3):
		VuL[i][0] = evYu2nt[i][1]
		VuL[i][1] = evYu2nt[i][2]
		VuL[i][2] = evYu2nt[i][0]
		VuR[i][0] = evYu2tn[i][1]
		VuR[i][1] = evYu2tn[i][2]
		VuR[i][2] = evYu2tn[i][0]

		VdL[i][0] = evYd2nt[i][1]
		VdL[i][1] = evYd2nt[i][2]
		VdL[i][2] = evYd2nt[i][0]
		VdR[i][0] = evYd2tn[i][1]
		VdR[i][1] = evYd2tn[i][2]
		VdR[i][2] = evYd2tn[i][0]

	YuNorm =((VuL.transpose()).dot(Yu)).dot(VuR)
	YdNorm =((VdL.transpose()).dot(Yd)).dot(VdR)
	print("#########################################")
	print("############Biunitarity##################")
	print("#########################################")
	print("YuNorm:")
	print(YuNorm)
	print("YdNorm:")
	print(YdNorm)
	print("*****************************************")
	
	return VuL, VuR, VdL, VdR

def compareMathPaper(): #function = 4
	VuLMath = biUnitaryDiag()[0]	
	VuRMath = biUnitaryDiag()[1]	
	VdLMath = biUnitaryDiag()[2]	
	VdRMath = biUnitaryDiag()[3]	

	VuLPaper = diagonality(1)[0]
	VuRPaper = diagonality(1)[2]
	VdLPaper = diagonality(0)[0]
	VdRPaper = diagonality(0)[2]
	
	print("Math",VdRMath)
	print("Paper",VdRPaper)
	print("Difference",abs(VdRMath)-abs(VdRPaper))


function = 0 #0sympy, 1diagonalize, 2unitarityOfV, 3biUnitarityDiag, 4compareMathPaper
upOrDown = 0 #0down, 1up
lOrR = 2 #0Vl, 2Vr
if(function == 0):
	sympy()
if(function == 1):
 	diagonality(upOrDown)
if(function == 2):
	unitarityOfV(upOrDown, lOrR)
if(function ==3):
	biUnitaryDiag()
if(function ==4):
	compareMathPaper()
