#!/usr/bin/python

import sys, os
import numpy as np
import re

if len(sys.argv) != 2:
	print "Usage: ./get_timings.py n_samples"
	exit(1)

n = int(sys.argv[1])

#run with 1 thread divergent
s = os.popen('./timing_test -n %d -t 100' % n).read()

a,b = s.split('---------');

l = []
for i in re.findall(r'\d+', a):
	l.append(float(i))
print "Div 1 thread mean:\t %f std: %f" % (np.mean(l), np.std(l))
with open('div.dat', 'w') as f:
	f.write(str(l)+'\n')

l = []
for i in re.findall(r'\d+', b):
	l.append(float(i))
print "No div 1 thread mean:\t %f std: %f" % (np.mean(l), np.std(l))
with open('nodiv.dat', 'w') as f:
	f.write(str(l)+'\n')

#run with no divergent threads
# l = []
# s = os.popen('./timing_test -n %d -t 100' % n).read()
# for i in re.findall(r'\d+\.\d+ ms', s):
# 	l.append(float(i.split()[0]))
# print "No div thread mean:\t %f std: %f" % (np.mean(l), np.std(l))
# with open('nodiv.dat', 'w') as f:
# 	f.write(str(l)+'\n')

