#!/usr/bin/python

import random, time

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
	#print("e:", e)
	#print(bits(e))
	
	for i in e_b[1:]:
		#print(i)
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
	return len(i_b), 2**len(i_b) #changes: more efficient

def REDC(R,N,N_,T,L): #changes: more efficient
	m = ((T & R) * N_) & R #changes: more efficient
	t = (T + m*N) >> L #changes: python 3 compatible
	if t >= N:
		t = t - N
		return t
	else:
		N = t - N
		return t

def CheckREDC(R,N,N_,T,L): #changes: more efficient
	m = ((T & R) * N_) & R #changes: more efficient
	t = (T + m*N) >> L #changes: python 3 compatible
	if t >= N:
		return 1
	else:
		return 0

def MontExp(mes,e,n):
	r = findR(n)[1] #changes: more efficient
	rmod = r - 1  #changes: more efficient	
	l = findR(n)[0] #changes: more efficient
	n_ = - modinv(n,r) & rmod #changes: more efficient
	r2 = (r << l) % n #changes: more efficient
	_a = REDC(rmod,n,n_,1*r2,l) #changes: more efficient
	_b = REDC(rmod,n,n_,mes*r2,l) #changes: more efficient
	for i in bits(e)[::-1]:
		if i == '1':
			_a = _a * _b
			_a = REDC(rmod,n,n_,_a,l) #changes: more efficient
		_b = _b * _b
		_b = REDC(rmod,n,n_,_b,l) #changes: more efficient
	_a = REDC(rmod,n,n_,_a,l) #changes: more efficient
	return _a

def MontSQMLadder(mes, e, n):
	r = findR(n)[1] #changes: more efficient
	rmod = r - 1  #changes: more efficient	
	l = findR(n)[0] #changes: more efficient
	n_ = - modinv(n,r) & rmod #changes: more efficient
	r2 = (r << l) % n #changes: more efficient
	_x1 = REDC(rmod,n,n_,mes*r2,l) #changes: more efficient
	_x2 = _x1 * _x1
	_x2 = REDC(rmod,n,n_,_x2,l)
	e_b = bits(e)
	for i in e_b[1:]:
		if i == '0':
			_x2 = _x1 * _x2
			_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
			_x1 = _x1 * _x1
			_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
		else:
			_x1 = _x1 * _x2
			_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
			_x2 = _x2 * _x2
			_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
	_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
	return _x1

def CheckDivExp(mes,e,n,bit):
	r = findR(n)[1] #changes: more efficient
	rmod = r - 1  #changes: more efficient	
	l = findR(n)[0] #changes: more efficient
	n_ = - modinv(n,r) & rmod #changes: more efficient
	r2 = (r << l) % n #changes: more efficient
	_x1 = REDC(rmod,n,n_,mes*r2,l) #changes: more efficient
	_x2 = _x1 * _x1
	_x2 = REDC(rmod,n,n_,_x2,l)
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
			d1_1 = CheckREDC(rmod,n,n_,_x2,l) #changes: more efficient
			_x1 = _x1 * _x1
			d1_2 = CheckREDC(rmod,n,n_,_x1,l) #changes: more efficient

			#simulate exp bit 1
			_x1 = _x1_temp
			_x2 = _x2_temp
			_x1 = _x1 * _x2
			d0_1 = CheckREDC(rmod,n,n_,_x1,l) #changes: more efficient
			_x2 = _x2 * _x2
			d0_2 = CheckREDC(rmod,n,n_,_x2,l) #changes: more efficient
			return ((d0_1, d0_2), (d1_1, d1_2))
		else:
			#print "[%d]=%s cont ..." % (c, i)
			if i == '0':
				_x2 = _x1 * _x2
				_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
				_x1 = _x1 * _x1
				_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
			else:
				_x1 = _x1 * _x2
				_x1 = REDC(rmod,n,n_,_x1,l) #changes: more efficient
				_x2 = _x2 * _x2
				_x2 = REDC(rmod,n,n_,_x2,l) #changes: more efficient
			c -= 1
	
def CalcDiv(bit0, bit1):
	di = lambda p0, p1: abs(p0[0]-p1[0]) + abs(p0[1]-p1[1])
	return list(map(di, bit0, bit1))

#find pairs that cause divergence when given exponent bit is 1
def FindDiv (num, mod, e, bit): 
	res = []
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] == 1 and CalcDiv(d1,d2)[1] == 0:
				res.append((r1,r2))
				print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
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
				print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break
			# else:
			# 	print "%d, %d, %s, %s" % (r1,r2,str(d1), str(d2))
	return res				

random.seed(time.time())
p = 32416189867
q = 32416189909
n = p*q 
phi = (p-1)*(q-1)
n_lambda = phi // egcd(p-1, q-1)[0] #changes: more efficient
e = 5
d = modinv(e, n_lambda) #changes: more efficient

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
R = findR(n)[1] #changes: more efficient
L = findR(n)[0] #changes: more efficient
R2 = (R << L) % n #changes: more efficient
a = 137
b = 262
c = a * b % n
am = toMont(a,R,n)
bm = toMont(b,R,n)
cm = toMont(c,R,n)
N_ = - modinv(n,R) % R
tmp1 = am * bm
tmp2 = REDC(R,n,N_,tmp1,L)
tmp3 = REDC(R,n,N_,tmp2,L)
print ("R: %d (%s) N_: %d a: %d am: %d b: %d bm: %d c: %d cm:%d" % (R, bits(R), N_, a, am, b, bm, c, cm))
print ("am*bm=%d, reduced: %d, double-reduced: %d " % (tmp1, tmp2, tmp3))

#Exponentation with Montgomery multiplication:
print ("")
c = MontExp(mes,e,n)
m2 = MontExp(c,d,n)
print ("MontSQM m: %d c: %d m_dec: %d" % (mes,c,m2))

#Exponentation with Montgomery multiplication:
print ("")
c = MontSQMLadder(mes,e,n)
m2 = MontSQMLadder(c,d,n)
print ("MontLadderSQM m: %d c: %d m_dec: %d" % (mes,c,m2))

print ("Secret exponent len: %d \nbits: %s" % (len(bits(d)), bits(d)))

d1 = CheckDivExp(1000,d,n,52)
d2 = CheckDivExp(1001,d,n,52)
print (d1,d2)

print (CalcDiv(d1,d2))
print (d)

print (FindDiv (10, n, d, 52))
print ("")
print (FindNoDiv (10, n, d, 52))


#print( "CUDA inputs: hex(n):%s, hex(N_):%s, hex(R):%s, hex(R2):%s, hex(RMOD):%s, bits(e):%s, bits(d):%s, L:%d" % (hex(n), hex(N_), hex(R), hex(R2), hex(R - 1), bits(e), bits(d), L))


hex_n = hex(n)[2:]
padding = 8 - (len(hex_n) % 8);
for i in range(padding):
	hex_n = "0" + hex_n;
	
hex_N_ = hex(n)[2:]
padding = 8 - (len(hex_n) % 8);
for i in range(padding):
	hex_N_ = "0" + hex_N_;
	
hex_R2 = hex(n)[2:]
padding = 8 - (len(hex_n) % 8);
for i in range(padding):
	hex_R2 = "0" + hex_R2;

print( "CUDA inputs: hex(n):%s, hex(N_):%s, hex(R2):%s, bits(e):%s, bits(d):%s, L:%d" % (hex_n, hex_N_, hex_R2, bits(e), bits(d), L))

