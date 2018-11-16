#!/usr/bin/python

import random, time
import numpy as np
import sys

def red(num,mod):
	return num % mod

def bits(i):
	return "{0:b}".format(i)

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

def sqm(mes,e,n):
	#right to left
	_a = 1
	_b = mes
	for i in bits(e)[::-1]:
		if i == '1':
			_a = _a * _b
			_a = red(_a, n)
		_b = _b * _b
		_b = red(_b, n)
	return _a

def SQMLadder(mes, e, n):
	x1 = mes 
	x2 = mes * mes
	x2 = red(x2,n)
	e_b = bits(e)
	for i in e_b[1:]:
		if i == '0':
			x2 = x1 * x2
			x2 = red(x2, n) 
			x1 = x1 * x1
			x1 = red(x1, n)
		else:
			x1 = x1 * x2
			x1 = red(x1, n)
			x2 = x2 * x2
			x2 = red(x2, n)
	return x1

def toMont(i, R, N):
	return (i * R) % N

def findR(i):
	i_b = "{0:b}".format(i)
	return 2**len(i_b)

def REDC(R,N,N_,T):
	m = ((T % R) * N_) % R
	t = (T + m*N) // R
	if t >= N:
		t = t - N
		return t
	else:
		N = t - N
		return t

def CheckREDC(R,N,N_,T):
	m = ((T % R) * N_) % R
	t = (T + m*N) // R
	if t >= N:
		return 1
	else:
		return 0

def MontExp(mes,e,n):
	r = findR(n)
	n_ = - modinv(n,r) % r
	_a = toMont(1,r,n)
	_b = toMont(mes,r,n)
	for i in bits(e)[::-1]:
		if i == '1':
			_a = _a * _b
			_a = REDC(r,n,n_,_a)
		_b = _b * _b
		_b = REDC(r,n,n_,_b)
	_a = REDC(r,n,n_,_a)
	return _a

def MontSQMLadder(mes, e, n):
	r = findR(n)
	n_ = - modinv(n,r) % r
	_x1 = toMont(mes,r,n)
	_x2 = _x1 * _x1
	_x2 = REDC(r,n,n_,_x2)
	e_b = bits(e)
	for i in e_b[1:]:
		if i == '0':
			_x2 = _x1 * _x2
			_x2 = REDC(r,n,n_,_x2)
			_x1 = _x1 * _x1
			_x1 = REDC(r,n,n_,_x1)
		else:
			_x1 = _x1 * _x2
			_x1 = REDC(r,n,n_,_x1)
			_x2 = _x2 * _x2
			_x2 = REDC(r,n,n_,_x2)
	_x1 = REDC(r,n,n_,_x1)
	return _x1

def CheckDivExp(mes,e,n,bit):
	r = findR(n)
	n_ = - modinv(n,r) % r
	_x1 = toMont(mes,r,n)
	_x2 = _x1 * _x1
	_x2 = REDC(r,n,n_,_x2)
	e_b = bits(e)
	if bit > len(e_b) - 2 :
		print ("Wrong bit!")
		exit(1)
	c = len(e_b) - 2
	for i in e_b[1:]:
		if bit == c:
			#print "Got ya! [%d]=%s" % (c, i)
			_x1_temp = _x1
			_x2_temp = _x2

			#simulate exp bit 0
			_x2 = _x1 * _x2
			d1_1 = CheckREDC(r,n,n_,_x2)
			_x1 = _x1 * _x1
			d1_2 = CheckREDC(r,n,n_,_x1)

			#simulate exp bit 1
			_x1 = _x1_temp
			_x2 = _x2_temp
			_x1 = _x1 * _x2
			d0_1 = CheckREDC(r,n,n_,_x1)
			_x2 = _x2 * _x2
			d0_2 = CheckREDC(r,n,n_,_x2)
			return ((d0_1, d0_2), (d1_1, d1_2))
		else:
			#print "[%d]=%s cont ..." % (c, i)
			if i == '0':
				_x2 = _x1 * _x2
				_x2 = REDC(r,n,n_,_x2)
				_x1 = _x1 * _x1
				_x1 = REDC(r,n,n_,_x1)
			else:
				_x1 = _x1 * _x2
				_x1 = REDC(r,n,n_,_x1)
				_x2 = _x2 * _x2
				_x2 = REDC(r,n,n_,_x2)
			c -= 1

def DivPattern(mes,e,n):
	pattern = []
	r = findR(n)
	n_ = - modinv(n,r) % r
	_x1 = toMont(mes,r,n)
	_x2 = _x1 * _x1
	_x2 = REDC(r,n,n_,_x2)
	e_b = bits(e)
	for i in e_b[1:]:
		_x1_temp = _x1
		_x2_temp = _x2

		if i == '0':

			#simulate exp bit 0
			__x2 = _x1 * _x2
			d1 = CheckREDC(r,n,n_,__x2)
			__x1 = _x1 * _x1
			d2 = CheckREDC(r,n,n_,__x1)
			
			_x2 = _x1 * _x2
			_x2 = REDC(r,n,n_,_x2)
			_x1 = _x1 * _x1
			_x1 = REDC(r,n,n_,_x1)
		else:
			__x1 = _x1 * _x2
			d1 = CheckREDC(r,n,n_,__x1)
			__x2 = _x2 * _x2
			d2 = CheckREDC(r,n,n_,__x2)
		
			_x1 = _x1 * _x2
			_x1 = REDC(r,n,n_,_x1)
			_x2 = _x2 * _x2
			_x2 = REDC(r,n,n_,_x2)

		pattern.append((d1,d2))

	_x1 = REDC(r,n,n_,_x1)
	#print _x1
	#print pattern
	return pattern
	
def CalcDiv(bit0, bit1):
	di = lambda p0, p1: abs(p0[0]-p1[0]) + abs(p0[1]-p1[1])
	return list(map(di, bit0, bit1))

#find pairs that cause divergence when given exponent bit is 0
def FindDiv (num, mod, e, bit): 
	res = []
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] == 1 and CalcDiv(d1,d2)[1] == 0:
				res.append((r1,r2))
				#print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break
			# else:
			# 	print "%d, %d, %s, %s" % (r1,r2,str(d1), str(d2))
	return res

#find pairs that cause NO divergence in given bit
def FindNoDiv (num, mod, e, bit): 
	res = []
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] == CalcDiv(d1,d2)[1] == 0:
				res.append((r1,r2))
				#print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break
			# else:
			# 	print "%d, %d, %s, %s" % (r1,r2,str(d1), str(d2))
	return res

def DivCoef(l, e, n):
	res = []
	for tup in l:
		c = 0
		p1 = DivPattern(tup[0], e, n)
		p2 = DivPattern(tup[1], e, n)
		#print p1

		for i, val in enumerate(p1):
			if val[0] != p2[i][0]:
				c += 1
			if val[1] != p2[i][1]:
				c += 1
		res.append(c)
	return res


random.seed(time.time())
p = 32416189867
q = 32416189909
n = p*q 
phi = (p-1)*(q-1)
n_lambda = phi // egcd(p-1, q-1)[0] #changes: more efficient
e = 5
d = modinv(e, phi)

e_b = bits(e)
d_b = bits(d)
n_b = bits(n)

print ("p: %d q: %d n: %d(%s) phi: %d e: %d(%s) d: %d(%s)" % (p,q,n,n_b,phi,e,e_b,d,d_b))

#encrypt:
mes = 12345
c = pow(mes, e, n)
#decrypt:
m1 = pow(c,d,n)
print ("Pow m: %d c: %d m_dec: %d" % (mes,c,m1))

#square-and-multiply
print ("")
c = sqm(mes,e,n)
m2 = sqm(c,d,n)
print ("SQM m: %d c: %d m_dec: %d" % (mes,c,m2))

#square-and-multiply with Montgomery ladder
print ("")
c = SQMLadder(mes,e,n)
m2 = SQMLadder(c,d,n)
print ("Ladder SQM m: %d c: %d m_dec: %d" % (mes,c,m2))

#Testing Montgomery multiplication:
R = findR(n)
a = 137
b = 262
c = a * b % n
am = toMont(a,R,n)
bm = toMont(b,R,n)
cm = toMont(c,R,n)
N_ = - modinv(n,R) % R
tmp1 = am * bm
tmp2 = REDC(R,n,N_,tmp1)
tmp3 = REDC(R,n,N_,tmp2)
print ("R: %d (%s) N_: %d a: %d am: %d b: %d bm: %d c: %d cm:%d" % (R, bits(R), N_, a, am, b, bm, c, cm))
print ("am*bm=%d, reduced: %d, double-reduced: %d " % (tmp1, tmp2, tmp3))

#Exponentation with Montgomery multiplication:
print ("")
c = MontExp(mes,e,n)
m2 = MontExp(c,d,n)
print ("MontSQM m: %d c: %d m_dec: %d" % (mes,c,m2))

#Exponentation with Montgomery multiplication:
print ("")
c = MontSQMLadder(mes,d,n)
m2 = MontSQMLadder(c,e,n)
print ("MontLadderSQM m: %d c: %d m_dec: %d" % (mes,c,m2))

print ("Secret exponent len: %d \nbits: %s" % (len(bits(d)), bits(d)))
print ("\n\n\n\n\n")

d1 = CheckDivExp(1000,d,n,52)
d2 = CheckDivExp(1001,d,n,52)
print (d1,d2)

print (CalcDiv(d1,d2))
print (d)

print(len(bits(d)))

print ("\n\n\n\n\n")
# print(DivPattern(746346390363684944272,d,n))
# print(enumerate(DivPattern(746346390363684944272,d,n)))

# print (FindDiv (10, n, d, 52))
# print ("")
# print (FindNoDiv (10, n, d, 52))



# l = FindDiv (1000, n, d, 52)
# print np.average(DivCoef(l, d, n))
# l = FindNoDiv(1000, n, d, 52)
# print np.average(DivCoef(l, d, n))

n_rep = 10
for i in range(2000,2001,100):
	probes = []
	for j in range(n_rep):
		l = FindDiv (i, n, d, 53)
		div = np.average(DivCoef(l, d, n))
		div_a =  np.average(div)
		l = FindNoDiv(i, n, d, 53)
		nodiv = np.average(DivCoef(l, d, n))
		nodiv_a = np.average(nodiv)
		delt = div_a - nodiv_a
		probes.append(delt)
		print ("Set size: %d div: %f nodiv: %f delta: %f" % (i, div_a, nodiv_a, delt))
	print ("Accuracy rate for size %d: %f%% avg:%f\n" % (i, float(sum(x > 0 for x in probes)) / n_rep * 100, np.average(probes)))
	sys.stdout.flush() 
print ("%f" % (float(sum(x > 0 for x in probes)) / n_rep * 100))

# calculate accuracy for sample size 500:
#for j in range(1000)