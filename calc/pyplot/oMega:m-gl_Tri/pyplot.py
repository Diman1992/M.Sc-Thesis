import numpy as np
import matplotlib.pyplot as plt

annrate = 4e-9 #1607.02457: 4e-9; 1512.01991: 7.5e-10; 1601.06590: 2.1e-8 (most convicing)
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

fig = plt.figure()
ax = plt.gca()
m = np.linspace(1,300,1000)

Mlf = 300
Mqf = 720

gofm = []
print("going to loop")
for i in m:
	gofm.append(anngToM(i,Mlf,Mqf))

	if(i%50<1):
		print(i)	


#ax.plot(m,tot)
#ax.plot(m,muon)
#ax.plot(m,top)
#ax.plot(m,WW)
#ax.plot(m,ann)

print("plot")

ax.plot(m,gofm)

ax.set_xlabel(r'$m_\chi$ / GeV')
ax.set_ylabel(r'$\sigma_{Ann}$ / pb')
ax.set_ylim(1,3)


ax.legend(loc='best')
plt.show()




