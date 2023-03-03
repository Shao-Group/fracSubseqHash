import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import MultipleLocator
import sys

stat = []
name = ['minimizer', 'alphabetical', 'pseudo_1_1', 'pseudo_11_1', 'pseudo_31_1']
label = ['Minimizer', 'Lexicographic', 'ABC, d=1', 'ABC, d=11', 'ABC, d=31']
m = ['s', '^', 'o', 'o', 'o', 'o', 'o', '*']

for i in name:
	f = open('./30/' + i)


	for p in range(15):
		f.readline()
		for j in range(20):
			s = f.readline().strip().split(',')

			p = 0.0
			if(len(s[0]) != 0):
				p = float(s[0])
			stat.append(p/100.0)		
		f.readline()

pr = []
f = open('./30/purerandom')

for p in range(6):
	f.readline()
	for j in range(10):
		s = f.readline().strip().split(',')

		pr.append(float(s[0])/100.0)
	f.readline()

titles = []
for i in range(15, 30):
	titles.append(str(i))

x = [1,2,3,4,5,6,7,8,9,10]
y = []

for z in range(15):
	plt.figure(figsize=(10,10))
	ax = plt.axes()

	for i in range(5):
		tmp = []
		for j in range(10):
			tmp.append(stat[i * 300 + z * 20 + j])

		#print(z,i,tmp)
		plt.plot(x, tmp, marker = m[i], label = label[i], linewidth=1.5)
		#plt.scatter(x, tmp, label = label[i])


	if(z <= 5):
		tmp = []
		for j in range(10):
			tmp.append(pr[z * 10 + j])
		plt.plot(x, tmp, marker = '*', label = 'fully random', linewidth=1.5)
		#plt.scatter(x, tmp)
	plt.grid()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(15) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(15) 
	ax.legend(fontsize=15)
	ax.set_xlabel('edit distance', fontsize=15)
	ax.set_ylabel('hash collision probability', fontsize=15)
	#ax.set_title()
	plt.savefig('./30/figures/k={0}.pdf'.format(titles[z]), bbox_inches='tight')
	plt.close()