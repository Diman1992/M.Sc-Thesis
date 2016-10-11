import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.interpolate import interp1d
fig = plt.figure()
ax1 = plt.subplot(111)


a = [] 
b = []
with open('out1.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		a.append(row[0])
		b.append(row[1])
ax1.plot(a,b,label=r'4,3')
for i in range(0,len(a)):
	a[i] = float(a[i])
	b[i] = float(b[i])

x = []
y = []
with open('out2.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		x.append(row[0])
		y.append(row[1])
ax1.plot(x,y,label=r'1,2')
for i in range(0,len(x)):
	x[i] = float(x[i])
	y[i] = float(y[i])


w = []
v = []
with open('out3.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		w.append(row[0])
		v.append(row[1])
ax1.plot(w,v,label=r'3,2')

c = []
d = []
with open('out4.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		c.append(row[0])
		d.append(row[1])
ax1.plot(c,d,label=r'3,2')

e = []
f = []
with open('out5.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		e.append(row[0])
		f.append(row[1])
ax1.plot(e,f,label=r'3,2')

g = []
h = []
with open('out6.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		g.append(row[0])
		h.append(row[1])
ax1.plot(g,h,label=r'3,2')


if(x[0]>x[-1]):
	np.fliplr(x)
	np.fliplr(y)
if(a[0]>a[-1]):
	np.fliplr(a)
	np.fliplr(b)


superxNew = np.linspace(min(a[0],x[0]),max(a[-1],x[-1]),len(x)+len(a))
superyNew = np.linspace(0,0,len(x)+len(a))



anew = np.linspace(a[0],a[-1],len(x))
xnew = np.linspace(x[0],x[-1],len(x))

bInt = interp1d(a,b)
yInt = interp1d(x,y)

bnew = bInt(anew)
#ax1.plot(anew,bnew,'')
ynew = yInt(xnew)
#ax1.plot(xnew,ynew,'o')
print(bnew, ynew)
ax1.fill_between(anew,ynew,bnew)

ax1.set_ylim(1,3)
ax1.grid('on')
ax1.legend(loc='best')

fig.savefig('g2m-gl_RepsLowHigh.pdf')



plt.show()

