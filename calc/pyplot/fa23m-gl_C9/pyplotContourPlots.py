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

ax1.plot(a,b,label=r'$C_9$high1')


x = []
y = []

with open('out2.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		x.append(row[0])
		y.append(row[1])

ax1.plot(x,y,label=r'$C_9$high2')
w = []
v = []


with open('out3.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		w.append(row[0])
		v.append(row[1])

ax1.plot(w,v,label=r'$C_9$low1')

ax1.set_ylim(1,3)
ax1.grid('on')
ax1.legend(loc='best')

fig.savefig('fa23m-gl_C9.pdf')




plt.show()
