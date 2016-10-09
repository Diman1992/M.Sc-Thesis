####
#Phenomenological functions and plots
#all entities with units in GeV
###

import numpy as np
import matplotlib.pyplot as plt

#processes = enumerate('bsprocess','g2mu','annihilation','directdetection')
#print( processes[0], processes[1])

weak = -4*1.16e-5/np.sqrt(2) * 0.22**2 * 1/137 / (4*np.pi)
cmToGeV2 = 2.57e27 #cm**2 = cmToGeV2 1/GeV**2
mmuon = 0.105
mW = 81
mH = 126
alphW = 1/137 * 1/0.23
Ncl = 1
Ncq = 3
fNuc = np.array([0.02,0.026,0.14,0,0,0])
fGlu = 1-np.sum(fNuc)
q2nd = np.array([0.22,0.11,0.026,0.019,0.012,0])
q2ndb = np.array([0.034,0.036,0.026,0.019,0.012,0])
qmass = np.array([2.3e-3, 4.8e-3, 95e-3, 1.275, 4.18, 173])
cQ = np.array([1,1,1,1.32,1.19,1])

print("fTG+q2 squared x 4/pi (xSec): ",((fGlu/9+(q2nd[4]+q2ndb[4])*3)/12)**2*4/np.pi)

#Constraints
C9low1f = -0.81
C9high1f = -0.51
C9low2f = -0.97
C9high2f = -0.37
bsmix = 2.5e-11
annrate = 4e-9 #1607.02457: 4e-9; 1512.01991: 7.5e-10; 1601.06590: 2.1e-8 (most convicing)
g2best = 2.87e-9
g2low = 2.07e-9
g2high = 3.57e-9

#model parameters
eps = 0.2
gqf = np.array([eps**3,eps**3,eps**2,eps**2,1,1])
mldiff = 200
mqdiff = 800
gq2f = gqf[2]
gq3f = gqf[4]
gl2f = 2.5

def tHad(m,M):
	return M**2/m**2

#g-2
def i(t):
	return 1/(12*(t-1)**4) * (2+3*t-6*t**2+t**3+6*t*np.log(t))
def ibar(t):
	return 1/t*i(1/t)
def preg2(g,m):
	return g**2*mmuon**2/(16*np.pi**2) * 1/(m+1)**2

def g2func(g,m,qf,qb):
	res = 0
	#if(len(qf) != len(qb)):
	#	print("qf and qb have different lengths")
	print("qf: ",np.sum(qf), " qb: ",np.sum(qb))
	
#	res += np.sum(qf)*i(tHad(m,m+mldiff)) + np.sum(qb)* 1/tHad(m,m+mldiff) * i(1/tHad(m,m+mldiff))
	res += np.sum(qf)*i(tHad(m,m+mldiff)) + np.sum(qb)* ibar(1/tHad(m,m+mldiff))
	return -preg2(g,m) * res

def gMutoMass(m,a,qf,qb):
	return np.sqrt(-a*16*np.pi**2 * m**2/mmuon**2 * (np.sum(qf)*i(tHad(m,m+mldiff)) + np.sum(qb)* ibar(1/tHad(m,m+mldiff)))**(-1))


qfB = np.array([0,-1])
qbB = np.array([-1,0])
qfA = np.array([0])
qbA = np.array([-1])


#Bsmixing and semileptonic
def Kshort(t):
	return (1-t+t**2*np.log(t))/(t-1)**2
def Gshort(t):
	return (1-t+t*np.log(t))/(t-1)**2
def Klong(tl,tq):
	return (Kshort(tl)-Kshort(tq))/(tl-tq)
def Glong(tl,tq):
	return (Gshort(tl)-Gshort(tq))/(tl-tq)
def Kder(t):
	return (-1+t+2*t*np.log(t))/(t-1)**2 - 2*(1-t+t**2*np.log(t))/(t-1)**3
def Gder(t):
	return (np.log(t))/(t-1)**2 - 2*(1-t+t*np.log(t))/(t-1)**3
def MuMuAddMaj(m,relPreG): #preK = 1, what is the relative prefactor of preG? Crivellin
	return Klong(tHad(m,m+mldiff),tHad(m,m+mqdiff)) + 1/2*relPreG*Glong(tHad(m,m+mldiff),tHad(m,m+mqdiff))
def MixAddMaj(m,relPreG): #preK = 1, what is the relative prefactor of preG? Crivellin
	return Kder(tHad(m,m+mqdiff)) + 1/2*relPreG*Gder(tHad(m,m+mqdiff))

def g23mumu(pMod,C,m,gl2,preG):
	return C*weak/pMod *m**2/gl2**2 * 1/MuMuAddMaj(m,preG)
def g23mix(pMod,m,preG):
	return bsmix*m**2/(pMod*MixAddMaj(m,preG))
def bla(m):
	return 0.04

preA = 1/(128*np.pi**2)
preB = 5/(384*np.pi**2)
preGA = 1
preGB = 1/5

#direct detection A
f1A = 1/16 * np.sum(fNuc * gqf**2) #tree channel to light quarks
f2A = 1/216 * np.sum(cQ*gqf**2)*fGlu #one loop with gluon (all quarks in loop included)
f3A = 3/16 * np.sum((q2nd+q2ndb)*gqf**2) #twist operators
def mRed(m,n):
	return m*n/(m+n)
def fdd(m,mqdiff):
	return (f1A+f2A+f3A)*m /(m+mqdiff)**4
def sigmaDDA(m,A):
	return 4/np.pi * mRed(m,A)**2 * (A * fdd(m,mqdiff))**2

#direct detection B
def b(x):
	return np.sqrt(1-x/4)
def ddatan(x):
	return np.arctan(2*b(x)/np.sqrt(x))
def ddatan2(x):
	return np.arctan(np.sqrt(4-x)/np.sqrt(x))
def gH(x):
	return -2/b(x)*(2+2*x-x**2)*ddatan(x) + 2*np.sqrt(x)*(2-x*np.log(x))
def gAV(x):
	return 1/(24*b(x)) * np.sqrt(x)*(8-x-x**2)*ddatan(x) -1/24 *x*(2-(3+x)*np.log(x))
def gT1(x):
	return 1/3*b(x) *(2+x**2)*ddatan(x) + 1/12 *np.sqrt(x) * (1-2*x-x*(2-x)*np.log(x)) #maybe its rather 1/b than b
def gT2(x):
	return 1/(4*b(x)) *(2-4*x+x**2)*ddatan(x) + 1/4 *np.sqrt(x) * (1-2*x-x*(2-x)*np.log(x)) #maybe its rather 1/b than b
def gB1(x):
	return -1/24 *np.sqrt(x) * (x*np.log(x)-2)+(x**2-2*x+4)/(24*b(x))*ddatan(x)
def gB31(x,y):
	return -x**(3/2)/(12*(y-x))  -x**(3/2)*y**2/(24*(y-x)**2)*np.log(y)  -x**(5/2)*(x-2*y)/(24*(y-x)**2)*np.log(x)  -x**(3/2)*np.sqrt(y)*(y+2)*np.sqrt(4-y)/(12*(y-x)**2)*ddatan2(y)  +x*(x**3-2*(y+1)*x**2+4*(y+1)*x+4*y)/(12*(y-x)**2*np.sqrt(4-x))*ddatan2(x)
def gB32(x,y):
	return -x**(3/2)*y/(12*(y-x)**2)  -x**(5/2)*y**2/(24*(y-x)**3)*np.log(y)  +x**(5/2)*y**2/(24*(y-x)**3)*np.log(x)  +x**(3/2)*np.sqrt(y)*(-6*y+x*y**2-2*x*y-2*x)/(12*(y-x)**3*np.sqrt(4-y))*ddatan2(y)  -x*y*(x**2*y-2*x*y-6*x-2*y)/(12*(y-x)**3*np.sqrt(4-x))*ddatan2(x)
def gW(x,y):
	return 2*gB1(x) + gB31(x,y) + cQ[4]*gB32(x,y)

def dq(m):
	return alphW**2/mW**2 * gAV(tHad(m,mW)) 
def f1B(m):
	return alphW**2/(4*mH**2*mW) * gH(tHad(m,mW)) * np.sum(fNuc)
def f2B(m):
	return 2/(9*mW)*alphW**2 * (1/mW**2 * gW(tHad(m,mW),tHad(m,qmass[5])) - 1/(12*mH**2) * gH(tHad(m,mW))) * fGlu
def f3B(m):
	return np.sum(q2nd+q2ndb) *alphW**2/mW**3 * (gT1(tHad(m,mW)) + gT2(tHad(m,mW)))
def sigmaDDsiB(m,A):
	return 4/np.pi * mRed(m,A)**2 * (A * (f1B(m)+f2B(m)+f3B(m)))**2
def sigmaDDsdB(m,A):
	return 12/np.pi * mRed(m,A)**2 * (dq(m)*(0.77-0.49-0.15))**2
print("f1B(fq): ",f1B(500), "; f2B(fG): ", f2B(500), "; f3B(g1+g2): ", f3B(500))
print("sigmaDDsdB: ",sigmaDDsdB(100,1))
print("sigmaDDsiB: ",sigmaDDsiB(100,1))
print("a/mw3 2 x 4/pi: ", (alphW**2/mW**3)**2*4/np.pi, " 2ndmom: ", np.sum(q2nd+q2ndb), " mW/18mH 2: ", mW**2/(18*mH**2))
chinow = 100
wnow = mW**2/chinow**2
tnow = qmass[5]**2/chinow**2
preG = 8*np.pi/9
preInG = alphW**2/(4*np.pi*mW)
pre = preG*preInG
print("pre: ", pre)
print("gW(mchi=100): ", gB1(wnow) + gB31(wnow,tnow), "gH(100): ",gH(wnow) )
print("fN in GeV-3 for mchi=100; 1:gW: ", pre*fGlu/mW**2*(2*gB1(wnow) + gB31(wnow,tnow)), " 2:gH: ", pre*fGlu*np.sum(cQ)/(3*mH**2)*gH(wnow), " 3:gTi: ",3/4*np.sum(q2nd+q2ndb)*alphW**2/mW**3*(gT1(wnow)+gT2(wnow)) , "4:fq: ", np.sum(fNuc)*alphW**2/(4*mH**2*mW)*gH(wnow))




fig = plt.figure()
ax = plt.gca()
m = np.linspace(qmass[5]/2+10,500,1000)
#m = np.linspace(50,500,1000)

"""
with open("outputGood.txt") as f:
	data = f.read()

data = data.split('\n')

x = [row.split(' ')[0] for row in data]
y = [row.split(' ')[1] for row in data]
ax.plot(x,y)
"""
"""
print(gMutoMass(g2best,100,0,-1))
ax.plot(m,gMutoMass(g2best,m,qfB,qbB))
ax.set_xlabel(r'$\mathcal{G}$')
"""

"""
#Direct Detection m-sigma
ax.plot(m,sigmaDDsiB(m,1)/cmToGeV2)
ax.set_xlabel(r'$m_\chi$ / GeV')
ax.set_ylabel(r'$\sigma_{SI}$ / cm$^2$')
"""

"""
#g-2 m-a
g2b_a = np.linspace(g2best,g2best,1000)
g2h_a = np.linspace(g2high,g2high,1000)
g2l_a = np.linspace(g2low,g2low,1000)
ax.plot(m,g2func(gl2f,m,qfA,qbA))
ax.plot(m,g2h_a,label=r'1\sigma')
ax.plot(m,g2b_a,label=r'best')
ax.plot(m,g2l_a)
ax.set_xlabel(r'$m_\chi$ / GeV')
ax.set_ylabel(r'$\Delta a_\mu$')
"""

"""
#BsPheno m-g2g3
#n = np.linspace(gq2f*gq3f,gq2f*gq3f,1000)
#ax.scatter(m,g2func(gl2f,m,qfB,qbB))
#ax.plot(m,(g23mumu(preA,C9low1f,m,gl2f,preGA)),label=r'$C_9 1\sigma$')
#ax.plot(m,(g23mumu(preA,C9high1f,m,gl2f,preGA)))
#ax.plot(m,(g23mumu(preA,C9low2f,m,gl2f,preGA)),label=r'$C_9 2\sigma$')
#ax.plot(m,(g23mumu(preA,C9high2f,m,gl2f,preGA)))
#ax.plot(m,(np.sqrt(g23mix(preB,m,preGA))),label=r'$C_{BB}$')
#ax.plot(m,n,label=r'\mathcal(G): $g^q_2 g^q_3$')
#ax.set_xlabel(r'$m_\chi$ / GeV')
#ax.set_ylabel(r'$g^q_2 g^q_3$')
"""
ax.legend(loc='best')
plt.show()


