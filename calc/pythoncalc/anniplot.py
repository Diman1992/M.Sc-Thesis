import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl

def crossSection(x,y,k,g):
	cC = cChi(x,y,k)
	cS = cSM(x,y)
	mX = mXi(x,y,k)
	G = gamma(x,y,k)
	if(g==1):
		mC = mChiGamma(x,y,k)
		return cC**2*cS**2/(8*np.pi*(4*mC**2 - mX**2)**2+G**2*mX**2*0) * 1e18
	else: 
		mC = mChi(x,y,k)
		return cC**2*cS**2/(8*np.pi*(4*mC**2 - mX**2)**2) * 1e18

def mChiGamma(x,y,k):
	mX = mXi(x,y,k)
	G = gamma(x,y,k)
	return np.sqrt((4*mX**6-np.sqrt(mX**6*(G**2*(1-16*mX**4)+mX**2)))/(16*mX**4-1)) 
	
def gamma(x,y,k):
	mX = mXi(x,y,k)
	cS = cSM(x,y)
	return mX*cS**2

def mChi(x,y,k):
	return np.sqrt(y)*k*x**2/np.sqrt(x**2+1.5) * 25.3

def mChibig(x,y,k):
	return np.sqrt(y)*k*x**2/np.sqrt(x**2-3.9) * 25.3

def mXi(x,y,k):
#	return 49
	return np.sqrt(y)*k*x*49

def cChi(x,y,k):
	mC = mChiGamma(x,y,k)
	mX = mXi(x,y,k)
	return mC/mX

def cSM(x,y):
	return y/lamb(x,y)

def lamb(x,y):
	return x*np.sqrt(y)*31e3

def relOmega(x,y,k):
	return y*mChi(x,y,k)**2/(mXi(x,y,k)*(4*mChi(x,y,k)**2-mXi(x,y,k)**2)) * 1.5e6 - lamb(x,y)

def func(x,n=1,M=1):
	res = x**2/((4*x**2*0-M**2)**2) * 1/(x*(x/(x+n))**2)
	return res

def Omega(x,y,k):
	return 0
	

fig1 = plt.figure()
ax = plt.gca()
#ax.set_xscale("log")
#ax.set_yscale("log")

h = np.linspace(1,7 ,num=700)
y = 1
k = np.array([0.5,1,2])
k = 1
per = 0.1*len(h)
print(4.7,gamma(4.7,y,k)**2*mXi(4.7,y,k)**2)


#for j in k:
for i in range(0,len(h)):
	if(i%per==0): print(i/per)
#		print(relOmega(h[i],y,k))
#		print(crossSection(h[i],y,k,1))
#		ax.scatter(h[i],mChiGamma(h[i],y,k),marker=".",s=2)
#		if(j==0.5):
#			ax.scatter(mChi(h[i],y,j),mXi(h[i],y,j),marker=".",s=1,color='red')
#			ax.scatter(mChiGamma(h[i],y,j),mXi(h[i],y,j),marker=".",s=2,color='red')
#		if(j==k[1]):
#			ax.scatter(mChi(h[i],y,j),mXi(h[i],y,j),marker="*",s=1,color='green')
#			ax.scatter(mChiGamma(h[i],y,j),mXi(h[i],y,j),marker="*",s=2, color='green')
#		if(j==k[2]):
#			ax.scatter(mChi(h[i],y,j),mXi(h[i],y,j),marker="o",s=1,color='blue')
#			ax.scatter(mChiGamma(h[i],y,j),mXi(h[i],y,j),marker="o",s=2,color='green')
#	ax.scatter(h[i],crossSection(h[i],y,k,1),marker=".",s=2)
#	ax.scatter(h[i],mXi(h[i],y,k))
#	ax.scatter(h[i],mChiGamma(h[i],y,k))
	ax.scatter(h[i],4*mChiGamma(h[i],y,k)**2-mXi(h[i],y,k)**2)

#		ax.scatter(h[i],func(h[i]), marker=".",s=2)
plt.show()
