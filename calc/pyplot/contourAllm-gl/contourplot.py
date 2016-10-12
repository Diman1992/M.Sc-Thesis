import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#Constraints
C9low1f = -0.71
C9high1f = -0.35
C9low2f = -0.91
C9high2f = -0.18
bsmix = 2.5e-11
annrate = 4e-9 #1607.02457: 4e-9; 1512.01991: 7.5e-10; 1601.06590: 2.1e-8 (most convicing)
g2best = 2.87e-9
g2low = 2.07e-9
g2high = 3.57e-9

weak = -4*1.16e-5/np.sqrt(2) * 0.22**2 * 1/137 / (4*np.pi)
cmToGeV2 = 2.57e27 #cm**2 = cmToGeV2 1/GeV**2
cmTopb2 = 10**20 #cm**2 = cmTopb2 * pb**2
mmuon = 0.105
mW = 81
mH = 126
alphE = 1/137
alphG = 0.1185 #at mZ
alphW = alphE /0.23
Ncl = 1
Ncq = 3
qmass = np.array([2.3e-3, 4.8e-3, 95e-3, 1.275, 4.18, 173])
preA = 1/(128*np.pi**2)
preB = 5/(384*np.pi**2)
preGrip = 7/(576*np.pi**2)
preQuint = 3/(265*np.pi**2)
preGA = 1
preGB = 1/5
preGG = 0
preGQ = 6/20
qfB = 1
qbB = 1

def tHad(m,M):
	return M**2/m**2
###########################
#Annihilation
###########################
def annPwave(m,M,mf):
	return 1/(12 *(M**2 - mf**2 + m**2)**2) * (8 *m**6 - 13 *m**4 * mf**2 - 5 *mf**2 *(M**2 - mf**2)**2 + 2 *m**2 *(4 *M**4 - 11 *M**2 *mf**2 + 5 *mf**4))
def annToFerm(m,M,g,mf,Nc):
	if (m<mf):
		return 0
	return Nc * g**4 * np.sqrt(1 - mf**2/m**2) /(32* np.pi *(M**2 - mf**2 + m**2)**2) *(mf**2 + annPwave(m, M, mf) *0.3**2)
def annToWW(m):
	if (m<mW):
		return 0
	return alphW**2/( 4*np.pi* m**2* (2 - mW**2/m**2)**2) *np.sqrt(1 - mW**2/m**2)**3
def annTotalB(m,Ml,Mq,gl):
	return annToFerm(m,Ml,gl,mmuon,Ncl) + annToFerm(m,Mq,1,qmass[5],Ncq) + annToWW(m) #+ 2*annToGG(m,Mq,1) + annToPhPh(m,Ml,gl,-1,Ncl)
def annTotalA(m,Ml,Mq,gl):
	return annToFerm(m,Ml,gl,mmuon,Ncl) + annToFerm(m,Mq,1,qmass[5],Ncq)
def anngToM(m,Ml,Mq):
	return ((annrate-annToFerm(m,Mq,1,qmass[5],Ncq)-annToWW(m))/annToFerm(m,Ml,1,mmuon,Ncl))**(1/4)

###########################
#Bsmixing and semileptonic
###########################
def Kshort(t):
	return (1-t+t**2*np.log(t))/(t-1)**2
def Gshort(t):
	return (1-t+t*np.log(t))/(t-1)**2
def Klong(tl,tq):
	return (Kshort(tl)-Kshort(tq))/(tl-tq)
def Glong(tl,tq):
	return (Gshort(tl)-Gshort(tq))/(tl-tq)
def MuMuAddMaj(m,Ml,Mq,relPreG): #preK = 1, what is the relative prefactor of preG? Crivellin
	return Klong(tHad(m,Ml),tHad(m,Mq)) + 2*relPreG*Glong(tHad(m,Ml),tHad(m,Mq))

def glmumu(pMod,C,m,Ml,Mq,preG):
	return (C*weak/pMod *m**2/(0.04) * 1/MuMuAddMaj(m,Ml,Mq,preG))**(1/2)
###########################
#Anomalous magnetic moment
###########################
def iloop(t):
	return 1/(12*(t-1)**4) * (2+3*t-6*t**2+t**3+6*t*np.log(t))
def iloopbar(t):
	return 1/t*iloop(1/t)
def preg2(g,m):
	return g**2*mmuon**2/(16*np.pi**2) * 1/(m)**2
def glda(m,a,Ml,qf,qb): #m against gl (therefore gl=1)
	res = preg2(1,m) * (qf*iloop(tHad(m,Ml)) + qbB*iloopbar(tHad(m,Ml)) )
	return np.sqrt(a / res)

print("qf*i + qb*ibar: ",qfB*iloop(tHad(100,300)) + qbB*iloopbar(tHad(100,300)) )
print("preg2(1,100): ",preg2(1,100))
print("glda(100,1,300,1,1): ",glda(100,1,300,1,1))
fig = plt.figure()
ax1 = plt.subplot(111)
m = np.linspace(1,400,1000)

Mlf = 300
Mqf = 720

plAnn = []
plBsmumuA = []
plBsmumuB = []
plBsmumuG = []
plBsmumuQ = []
plDaHigh = []
plDaLow = []
print("going to loop")
for i in m:
	plAnn.append(anngToM(i,Mlf,Mqf))
	plBsmumuA.append(glmumu(preA,C9high1f,i,Mlf,Mqf,preGA))
	plBsmumuB.append(glmumu(preB,C9high1f,i,Mlf,Mqf,preGB))
	plBsmumuG.append(glmumu(preGrip,C9high1f,i,Mlf,Mqf,preGG))
	plBsmumuQ.append(glmumu(preQuint,C9high1f,i,Mlf,Mqf,preGQ))
	plDaHigh.append(glda(i,g2high,Mlf,qfB,qbB))
	plDaLow.append(glda(i,g2low,Mlf,qfB,qbB))
	if(i%50<1):
		print(i)	

print("plot")

#ax.plot(m,plAnn)
ax1.plot(m,plBsmumuA,label=r'Singlet')
ax1.plot(m,plBsmumuB,label=r'Triplet')
ax1.plot(m,plBsmumuG,label=r'Quadruplet')
ax1.plot(m,plBsmumuQ,label=r'Quintuplet')
#ax.plot(m,plDaHigh)
#ax.plot(m,plDaLow)
matplotlib.rcParams.update({'font.size': 15})

ax1.set_xlabel(r'$m_\chi$ / GeV')
ax1.set_ylabel(r'$g_2^l$')

ax1.grid('on')
ax1.set_ylim(1.0,3)



ax1.legend(loc='best')
fig.savefig('BsmumuReps.pdf')
plt.show()






