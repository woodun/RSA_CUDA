#!/usr/bin/python2.7

import matplotlib
#matplotlib.use('Agg')

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

Z_CUT = 1000

def getf(s):                  
	return [float(x) for x in re.findall(r'\d+\.\d+', s)]

if len(sys.argv) == 1:
	sys.argv.append('div.dat')
	sys.argv.append('nodiv.dat')

with open(sys.argv[1]) as f:
	a = f.read()

with open(sys.argv[2]) as f:
	b = f.read()

a = getf(a)
b = getf(b)



#filter
z_avg = np.mean(stats.zscore(a))
l_before = len(a)
a = zip(a, stats.zscore(a))
a = filter(lambda x: False if np.abs(x[1])>= Z_CUT else True, a)
print "In %s z-score average %f filtered out %d samples" % (sys.argv[1], z_avg, l_before - len(a))

z_avg = np.mean(stats.zscore(b))
l_before = len(b)
b = zip(b, stats.zscore(b))
b = filter(lambda x: False if np.abs(x[1])>= Z_CUT else True, b)
print "In %s z-score average %f filtered out %d samples" % (sys.argv[2], z_avg, l_before - len(b))

a = zip(*a)[0]
b = zip(*b)[0]

#print a, len(a)
#print b, len(b)

print "%s mean: %f std: %f" % (sys.argv[1], np.mean(a), np.std(a))
print "%s mean: %f std: %f" % (sys.argv[2], np.mean(b), np.std(b))
print stats.ttest_ind(a,b)

import numpy as np
import matplotlib.pyplot as plt

fig, ax0 = plt.subplots()

ax0.hist(a, 10, range=[min(a), max(a)], facecolor='r', alpha=0.5, label=sys.argv[1])
ax0.hist(b, 10, range=[min(b), max(b)], facecolor='b', alpha=0.5, label=sys.argv[2])
ax0.legend()

# t1 = list(np.arange(0, 0+len(b), 1))
# t2 = list(np.arange(1, 1+len(b), 1))

# ax0.plot(t2, a, color='r', alpha=0.25, marker=',', linestyle='None', markersize='0.5', label=sys.argv[1])
# ax0.plot(t1, b, color='b', alpha=0.25, marker=',', linestyle='None', markersize='0.5', label=sys.argv[2])


# axes = plt.gca()
# axes.set_ylim([1.50,2.5])
# ax0.legend()

#fig.tight_layout()
fig.savefig('out.png')
plt.show()