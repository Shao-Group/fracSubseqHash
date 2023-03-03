import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import MultipleLocator
import sys

stat = []
titles = []
name = ['ABt0.txt', 'ABt785.txt', 'ABt-k-k2n.txt']
label = ['ABC, d=1', 'AB, t=0.785', 'AB, t=1-k/n']
m = ['s', '^', 'o', 'o', 'o', 'o', 'o', '*']
for i in name:
	f = open('./n20/' + i)


	for p in range(4):
		titles.append(f.readline().strip())
		for j in range(10):
			s = f.readline().strip().split(',')
			p = 0.0
			if(len(s[0]) != 0):
				p = float(s[0])
			stat.append(p/100.0)
		f.readline()


x = [1,2,3,4,5,6,7,8,9,10]
y = []

for z in range(4):
	plt.figure(figsize=(10,10))
	ax = plt.axes()

	for i in range(2):
		tmp = []
		for j in range(10):
			tmp.append(stat[i * 40 + z * 10 + j])

		plt.plot(x, tmp, marker = m[i], label = label[i], linewidth=1.5)
		#plt.scatter(x, tmp, label = label[i])


	plt.grid()

	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(15) 
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(15) 

	ax.legend(fontsize=15)
	ax.set_xlabel('edit distance', fontsize=15)
	ax.set_ylabel('hash collision probability', fontsize=15)
	#ax.set_title('N = 20, k = ' + titles[z])
	plt.savefig('./n20/figures/k={0}.pdf'.format(titles[z]), bbox_inches='tight')
	plt.close()
