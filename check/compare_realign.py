"""Compare the realigner outputs for different
version
ZZJ
2019.5.27
"""

import sys
import pysam



def read_as_score(bfile):
	s1 = {}
	with pysam.Samfile(bfile, 'rb') as bam:
		i = 0
		for r1 in bam:
			i += 1
			#if i>30000:
			#	break
			s1[(r1.qname, r1.rname, r1.pos)] = r1.opt('AS')
	return s1


def plot_scatter(new, old):
	import matplotlib.pyplot as plt
	import seaborn as sns
	ax = sns.jointplot(new, old, kind="reg")
	ax.set_axis_labels('New', 'Old')
	plt.savefig('realign_check.png')

def compare():
	s1 = read_as_score('new_out/realigned.sorted.bam')
	s2 = read_as_score('old_out/realigned.sorted.bam')
	k = list([x for x in s1 if x in s2])
	old = []
	new = []
	print("ID\tnew\told\n")
	for k_ in k:
		print("%s\t%s\t%s\n"%(k_, s1[k_], s2[k_] ) )
		new.append(s1[k_])
		old.append(s2[k_])
	plot_scatter(new, old)


if __name__ == '__main__':
	compare()