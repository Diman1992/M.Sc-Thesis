import numpy as np
import matplotlib.pyplot as plt

mmuon = 0.105
Ncl = 1
g2best = 2.87e-9
g2low = 2.07e-9
g2high = 3.57e-9

gl2f = 2.5
Mlf = 500

def tHad(m,M):
	return M**2/m**2
def i(t):
	return 1/(12*(t-1)**4) * (2+3*t-6*t**2+t**3+6*t*np.log(t))
def ibar(t):
	return 1/t*i(1/t)
def preg2(g,m):
	return g**2*mmuon**2/(16*np.pi**2) * 1/(m+1)**2
def g2func(g,m,M,qf,qb):
	res = np.sum(qf)*i(tHad(m,M)) + np.sum(qb)* ibar(1/tHad(m,M))
	return -preg2(g,m) * res

qfB = np.array([0,-1])
qbB = np.array([-1,0])
qfA = np.array([0])
qbA = np.array([-1])

fig = plt.figure()
ax = plt.gca()
m = np.linspace(0,400,1000)

g2b_a = np.linspace(g2best,g2best,1000)
g2h_a = np.linspace(g2high,g2high,1000)
g2l_a = np.linspace(g2low,g2low,1000)
ax.plot(m,g2func(gl2f,m,Mlf,qfA,qbA))
ax.plot(m,g2h_a)
ax.plot(m,g2b_a,label=r'best')
ax.plot(m,g2l_a)
ax.set_ylim(0,g2best*2)
ax.set_xlabel(r'$m_\chi$ in GeV')
ax.set_ylabel(r'$\Delta a_\mu$')

ax.legend(loc='best')
plt.show()


