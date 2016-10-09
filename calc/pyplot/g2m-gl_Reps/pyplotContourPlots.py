import matplotlib.pyplot as plt
import numpy as np
import csv
fig = plt.figure()
ax1 = plt.subplot(111)


a = []
b = []


with open('output3.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		a.append(row[0])
		b.append(row[1])

ax1.plot(a,b,label=r'4,3')


x = []
y = []

with open('output1.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		x.append(row[0])
		y.append(row[1])

ax1.plot(x,y,label=r'1,2')
w = []
v = []


with open('output2.txt','r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		w.append(row[0])
		v.append(row[1])

ax1.plot(w,v,label=r'3,2')

ax1.set_ylim(1,3)
ax1.grid('on')
ax1.legend(loc='best')


fig.savefig('g2m-gl_Reps.pdf')



plt.show()
