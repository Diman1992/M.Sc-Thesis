import numpy as np
import matplotlib.pylab as plt

def onlyL(t):
	return t**2*(2*((t-3)/(4*(t-1)**2)+(np.log(t))/(2*(t-1)**3)))

def onlyR(t):
	return t**2*((6*t-2)/(4*(t-1)**2) + (2*np.log(t))/(2*(t-1)**3) + (-1)/(t-1) + (-4*t*np.log(t))/(4*(t-1)**2))

def both(t):
	return onlyL(t) + onlyR(t)
	
def preconst(f,mB,lmu):
	return lmu**2/(16*np.pi**2)*f

fig1 = plt.figure()
ax=plt.gca()
h = np.linspace(0.01,0.9,100)
lmu = 1e-2
for i in range(0,len(h)):
	mB = 100/h[i]
#	ax.scatter(h[i],preconst(onlyL(h[i]),mB,lmu))
	print(h[i],preconst(onlyL(h[i]),mB,lmu))
#plt.show()
