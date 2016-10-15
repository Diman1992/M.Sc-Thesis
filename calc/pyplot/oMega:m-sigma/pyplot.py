import numpy as np
import matplotlib.pyplot as plt

annrate = 7.7e-9 #1607.02475: 7.7e-9; 1512.01991: 7.5e-10; 1601.06590: 2.1e-8 (most convicing)
weak = -4*1.16e-5/np.sqrt(2) * 0.22**2 * 1/137 / (4*np.pi)
cm2ToGeV2 = 2.57e27 #cm**2 = cmToGeV2 1/GeV**2
cmTopb = 10**20 #cm**2 = cmTopb2 * pb**2
cm3s1ToGeV2 = 7.7e17
mmuon = 0.105
mW = 81
mH = 126
alphE = 1/137
alphG = 0.1185 #at mZ
alphW = alphE /0.23
Ncl = 1
Ncq = 3
qmass = np.array([2.3e-3, 4.8e-3, 95e-3, 1.275, 4.18, 173])

Mlf = 300
Mqf = 720

def Li(z,n):
	res = 0
	for i in range(1,n):
		res+=z**i/i**2
	return res



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

def annToGG(m,Mq,g):
	n = 100
	return 2*alphG**2*g**4/(256*np.pi**3 *m**2) * (Li(-m**2/Mq**2,n) - Li(m**2/Mq**2,n))**2

def annToPhPh(m,M,g,Qf,Nc):
	n = 100
	return Nc**2 *Qf**4 * alphE**2 *g**4 /(256*np.pi**3*m**2) * (Li(-m**2/M**2,n) - Li(m**2/M**2,n))**2

def annTotalB(m,Ml,Mq,gl):
	return annToFerm(m,Ml,gl,mmuon,Ncl) + annToFerm(m,Mq,1,qmass[5],Ncq) + annToWW(m) #+ 2*annToGG(m,Mq,1) + annToPhPh(m,Ml,gl,-1,Ncl)

def annTotalA(m,Ml,Mq,gl):
	return annToFerm(m,Ml,gl,mmuon,Ncl) + annToFerm(m,Mq,1,qmass[5],Ncq)

gl2f = 1.8

fig = plt.figure()
ax = plt.gca()
m = np.linspace(0,250,1000)
tot = []
muon = []
top = []
WW = []
GG = []
PhPh = []
ann = []
for i in m:
	tot.append(annTotalB(i,Mlf,Mqf,gl2f)/cmTopb*cm2ToGeV2 )
	muon.append(annToFerm(i,Mlf,gl2f,mmuon,Ncl)/cmTopb*cm2ToGeV2 )
	top.append(annToFerm(i,Mqf,1,qmass[5],Ncq)/cmTopb*cm2ToGeV2)
	WW.append(annToWW(i)/cmTopb*cm2ToGeV2)
	ann.append(annrate/cmTopb*cm2ToGeV2)
#	GG.append(annToGG(i,Mqf,1))
#	PhPh.append(annToPhPh(i,Mlf,2.1,-1,Ncl))
#	if(i%10<1):
#		print(i,annToWW(i))	

ax.plot(m,tot,label=r'altogether')
ax.plot(m,muon,label=r'$\bar \mu\mu$')
ax.plot(m,top,label=r'$\bar t t$')
ax.plot(m,WW,label=r'$W^+W^-$')
ax.plot(m,ann,label=r'thermal relic')
#ax.plot(m,GG)
#ax.plot(m,PhPh)

ax.set_xlabel(r'$m_\chi$ / GeV')
ax.set_ylabel(r'$\sigma_{Ann}$ / pb')

ax.grid('on')

ax.legend(loc='best')
fig.savefig('annTriplet1.8.450900.pdf')
plt.show()




