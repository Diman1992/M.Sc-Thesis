import matplotlib.pyplot as plt
import numpy as np
import csv
fig = plt.figure()
ax1 = plt.subplot(111)

a = []
b = []
with open('out1.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		a.append(row[0])
		b.append(row[1])

ax1.plot(a,b,label=r'Singlet')

x = []
y = []
with open('out2.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		x.append(row[0])
		y.append(row[1])

ax1.plot(x,y,label=r'Triplet')

w = []
v = []
with open('out3.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		w.append(row[0])
		v.append(row[1])

ax1.plot(w,v,label=r'Quadruplet')

c = []
d = []
with open('out4.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		c.append(row[0])
		d.append(row[1])

ax1.plot(c,d,label=r'Quintuplet')

ax1.set_xlabel(r'$m_\chi$ / GeV')
ax1.set_ylabel(r'$g_2^l$')
plt.title('$g_2^l(m_\chi)$ for different $\chi$ irreps')

ax1.set_ylim(1.8,3)
ax1.grid('on')
ax1.legend(loc='best')

fig.savefig('fa23m-gl_Reps.pdf')
plt.show()
