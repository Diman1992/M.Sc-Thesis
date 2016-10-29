import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib import collections  as mc

#Constraints
C9low1f = -0.71
C9high1f = -0.35
C9low2f = -0.91
C9high2f = -0.18

C9low1f = -0.81
C9high1f = -0.51
C9low2f = -0.97
C9high2f = -0.37

bsmix = 0.34e-11
annrate = 7.7e-9 #1607.02457: 4e-9; 1512.01991: 7.5e-10; 1601.06590: 2.1e-8 (most convicing)
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
pre34 = 2/3*preA
preGA = 1
preGB = 1/5
preGG = 0
preGQ = 6/20
preG34 = 1/3
qfB = 1
qbB = 1
qf34 = 0
qb34 = 3
bsA = 0
bsB = 120
bs34 = 30 #bs34=bsB for Mlq=770

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
def anngToMB(m,Ml,Mq):
	return ((annrate-annToFerm(m,Mq,1,qmass[5],Ncq)-annToWW(m))/annToFerm(m,Ml,1,mmuon,Ncl))**(1/4)
def anngToMA(m,Ml,Mq):
	return ((annrate-annToFerm(m,Mq,1,qmass[5],Ncq))/annToFerm(m,Ml,1,mmuon,Ncl))**(1/4)


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
def Kder(t):
	return (-1+t+2*t*np.log(t))/(t-1)**2 - 2*(1-t+t**2*np.log(t))/(t-1)**3
def Gder(t):
	return (np.log(t))/(t-1)**2 - 2*(1-t+t*np.log(t))/(t-1)**3
def MuMuAddMaj(m,Ml,Mq,relPreG): #preK = 1, what is the relative prefactor of preG? Crivellin
	return Klong(tHad(m,Ml),tHad(m,Mq)) + 1*relPreG*Glong(tHad(m,Ml),tHad(m,Mq))
def MixAddMaj(m,relPreG): #preK = 1, what is the relative prefactor of preG? Crivellin
	return Kder(tHad(m,Mqf)) + 1*relPreG*Gder(tHad(m,Mqf))

def g23mix(pMod,m,preG):
	return bsmix*m**2/(pMod*MixAddMaj(m,preG))
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
	res = preg2(1,m) * (qf*iloop(tHad(m,Ml)) + qb*iloopbar(tHad(m,Ml)) )
	return np.sqrt(a / res)

print("qf*i + qb*ibar: ",qfB*iloop(tHad(100,300)) + qbB*iloopbar(tHad(100,300)) )
print("preg2(1,100): ",preg2(1,100))
print("glda(100,1,300,1,1): ",glda(100,1,300,1,1))
fig = plt.figure()
ax1 = plt.subplot(111)

matplotlib.rcParams.update({'font.size': 15})
ax1.set_xlabel(r'$m_\chi$ / GeV')
ax1.set_ylabel(r'$g_2^l$')
#ax1.grid('on')



Mlf = 300
Mqf = 720	

def contourAll(j,bs_,preM_,preG_,qf_,qb_,annx1_,annx2_,mmax_,ymin_,ymax_): #Reps 342
	m = np.linspace(1,mmax_,1000)
	ann_a = np.linspace(annx1_,annx2_,100)
	bsB = bs_
	pre = preM_
	pG = preG_
	qf = qf_
	qb = qb_
	plAnn = []
	plAnn2 = []
	plBsmumulow2 = []
	plBsmumuhigh2 = []
	plBsmumulow = []
	plBsmumu = []
	plDaHigh = []
	plDaLow = []
	print("going to loop")
	for i in ann_a: #for ann segment is solution
		plAnn2.append(anngToMB(i,Mlf,Mqf))
	for i in m:
		if(j==122):
			plAnn.append(anngToMA(i,Mlf,Mqf))
		else:
			plAnn.append(anngToMB(i,Mlf,Mqf))
		plBsmumuhigh2.append(glmumu(pre,C9high2f,i,Mlf,Mqf,pG))
		plBsmumulow2.append(glmumu(pre,C9low2f,i,Mlf,Mqf,pG))
		plBsmumulow.append(glmumu(pre,C9low1f,i,Mlf,Mqf,pG))
		plBsmumu.append(glmumu(pre,C9high1f,i,Mlf,Mqf,pG))
		plDaHigh.append(glda(i,g2high,Mlf,qf,qb))
		plDaLow.append(glda(i,g2low,Mlf,qf,qb))
		if(i%50<1):
			print(i)
		

	print("plot")
	test = plt.bar(10,20,width=20,color='b',bottom=50,alpha=0.5)
	extra = Rectangle((10, 2), 10, 10, fc="r", fill=True, edgecolor='none', linewidth=0)
	ax1.plot(m,plAnn,linestyle='dotted',label=r'$\Omega_{{DM}}$')
	ax1.plot(ann_a,plAnn2,color='b')
	#ax1.plot(m,plBsmumu)#,label=r'$C_9^{1}(1\sigma)$')
	#ax1.plot(m,plBsmumulow)
	#ax1.plot(m,plBsmumuhigh2)#,label=r'$C_9^{{1}}(2\sigma)$')
	#ax1.plot(m,plBsmumulow2)
	#ax1.plot(m,plDaHigh)#,label=r'high')
	#ax1.plot(m,plDaLow)#,label=r'low')
	#plt.plot((bsA,bsA),(1,3),'k-',)

##################################
# 122 patches: 1sigma: c:0050FF, a:0.6; 2sigma: c:B8CC8F, a:0.5; da: c:ff0000, a:0.6; bs: c:yellow, a:0.2
#     filling: 1sigma: c:0050FF, a:0.5; 2sigma: c:afaaaa, a:0.5; da: c:ff0000, a:0.6; bs: c:yellow, a:0.1
# 322 patches: 1sigma: c:ccfeaf, a:0.6; 2sigma: c:afaaaa, a:0.5; da: c:ff0000, a:0.6; bs: c:yellow, a:0.2
#     filling: 1sigma: c:0050FF, a:0.5; 2sigma: c:afaaaa, a:0.5; da: c:ff0000, a:0.6; bs: c:yellow, a:0.2
# 342 patches: 1sigma: c:0050FF, a:0.6; 2sigma: c:afaaaa, a:0.5; da: c:ff0000, a:0.6; bs: c:yellow, a:0.2
#     filling: 1sigma: c:0050FF, a:0.5; 2sigma: c:afaaaa, a:0.5; da: c:ff0000, a:0.6; bs: c:yellow, a:0.1
##################################
	ax1.add_patch(patches.Rectangle((0, 0),0.1,0.1, label=r'$C_9(1\sigma)$', facecolor = '#0050FF', alpha = 0.6	)) 
	ax1.add_patch(patches.Rectangle((0, 0),0.1,0.1, label=r'$C_9(2\sigma)$', facecolor = "#afaaaa", alpha = 0.5	))
	if(j==342):
		ax1.add_patch(patches.Rectangle((0, 0),0.1,0.1, label=r'$\Delta a_\mu$', facecolor = "#ff0000", alpha = 0.6	))
	ax1.add_patch(patches.Rectangle((0, 0),0.1,0.1, label=r'$\Delta m_s$', facecolor = "yellow", alpha = 0.2	))
	ax1.axvspan(bsB, mmax_, alpha=0.2, color='yellow')
	plt.fill_between(m, plBsmumuhigh2, plBsmumulow2, color='#afaaaa', alpha='0.5')
	plt.fill_between(m, plBsmumu, plBsmumulow, color='#0050FF', alpha='0.5')
	plt.fill_between(m, plDaHigh, plDaLow, color='#ff0000', alpha='0.6')
	ax1.set_ylim(ymin_,ymax_)
	#ax1.legend(loc=0)	
	ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.savefig('contour.pdf')

def bsmumuReps():
	m = np.linspace(1,250,1000)
	plBsmumuA = []
	plBsmumuB = []
	plBsmumuG = []
	plBsmumu34 = []
	print("going to loop")
	for i in m:
		plBsmumuA.append(glmumu(preA,C9high1f,i,Mlf,Mqf,preGA))
		plBsmumuB.append(glmumu(preB,C9high1f,i,Mlf,Mqf,preGB))
		plBsmumuG.append(glmumu(preGrip,C9high1f,i,Mlf,Mqf,preGG))
		plBsmumu34.append(glmumu(pre34,C9high1f,i,Mlf,Mqf,preG34))
		if(i%50<1):
			print(i)
		

	ax1.plot(m,plBsmumu34,label=r'$C_9^{342}(1\sigma)$')
	ax1.plot(m,plBsmumuA,linestyle='', marker='o',markevery=10,markersize=1,label=r'$C_9^{122}(1\sigma)$')
	ax1.plot(m,plBsmumuG,linestyle='-.',label=r'$C_9^{433}(1\sigma)$')
	ax1.plot(m,plBsmumuB,linestyle=':',label=r'$C_9^{322}(1\sigma)$')
	ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax1.set_ylim(1.4,3)
	
def bsmixReps():
	ax1.set_ylabel(r'$g_2^{q*}g_3^q$')
	m = np.linspace(1,200,1000)
	plBsmixA = []
	plBsmixAM = []
	plBsmixB = []
	plBsmixBM = []
	g2g3 = []
	print("going to loop")
	for i in m:
		plBsmixAM.append(np.sqrt(g23mix(preA,i,preGA)))
		plBsmixBM.append(np.sqrt(g23mix(preB,i,preGB)))
		plBsmixA.append(np.sqrt(g23mix(preA,i,0)))
		plBsmixB.append(np.sqrt(g23mix(preB,i,0)))
		g2g3.append(0.04)
		if(i%50<1):
			print(i)
	#ax1.add_patch(patches.Rectangle((0, 0),0.1,0.1, label=r'$\Delta m_s^1$', facecolor = '#00ff88', alpha = 0.5	)) 
	#ax1.add_patch(patches.Rectangle((0, 0),0.1,0.1, label=r'$\Delta m_s^3$', facecolor = "#ff0088", alpha = 0.5	))
	#plt.fill_between(m, plBsmixA, plBsmixAM, color='#00ff88', alpha='0.5')
	#plt.fill_between(m, plBsmixB, plBsmixBM, color='#ff0088', alpha='0.5')
	ax1.plot(m,plBsmixAM,linestyle='-.',label=r'$\Delta m_s^{1M}$')
	ax1.plot(m,plBsmixA,linestyle='-.',label=r'$\Delta m_s^{1D}$')
	ax1.plot(m,plBsmixBM,linestyle=':',label=r'$\Delta m_s^{3M}$')
	ax1.plot(m,plBsmixB,linestyle=':',label=r'$\Delta m_s^{3D}$')
	ax1.plot(m,g2g3,label=r'$\mathcal{G}: g^{q*}_2 g^q_3$')
	ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax1.set_ylim(0.035,0.06)

j = 1342
i = 0

if(j==122):
	contourAll(j,bsA,preA,preGA,0,1,53.2,81.2,150,1.5,3)
if(j==322):
	contourAll(j,bsB,preB,preGB,1,1,bsB,148,250,1.4,2.4)
if(j==342): #Mlq = 770
	contourAll(j,bs34,pre34,preG34,0,3,30,52.5,150,2,3)

#j,bs_,preM_,preG_,qf_,qb_,annx1_,annx2_,mmax_,ymin_,ymax_


if(i == 0): #Bsmumu diff reps
	bsmumuReps()
if(i == 1): #bsmix diff reps
	bsmixReps()


#ax.plot(m,plDaHigh)
#ax.plot(m,plDaLow)

 

plt.show()






