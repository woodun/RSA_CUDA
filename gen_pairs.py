#!/usr/bin/python

import random, time, sys

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

def findR(i):
	i_b = "{0:b}".format(i)
	return len(i_b), 2**len(i_b) 

def REDC(R,N,N_,T,L): 
	m = ((T & R) * N_) & R 
	t = (T + m*N) >> L 
	if t >= N:
		t = t - N
		return t
	else:
		N = t - N
		return t

def CheckREDC(R,N,N_,T,L): 
	m = ((T & R) * N_) & R 
	t = (T + m*N) >> L
	if t >= N:
		return 1
	else:
		return 0

def CheckDivExp(mes,e,n,bit):
	r = findR(n)[1] 
	rmod = r - 1  	
	l = findR(n)[0] 
	n_ = - modinv(n,r) & rmod 
	r2 = (r << l) % n 
	_x1 = REDC(rmod,n,n_,mes*r2,l) 
	_x2 = _x1 * _x1
	_x2 = REDC(rmod,n,n_,_x2,l)
	e_b = bits(e)
	if bit > len(e_b) - 2 :
		print ("Wrong bit!")
		exit(1)
	c = len(e_b) - 2
	for i in e_b[1:]:
		if bit == c:			
			_x1_temp = _x1
			_x2_temp = _x2

			#simulate exp bit 0
			_x2 = _x1 * _x2
			d1_1 = CheckREDC(rmod,n,n_,_x2,l) 
			_x1 = _x1 * _x1
			d1_2 = CheckREDC(rmod,n,n_,_x1,l) 

			#simulate exp bit 1
			_x1 = _x1_temp
			_x2 = _x2_temp
			_x1 = _x1 * _x2
			d0_1 = CheckREDC(rmod,n,n_,_x1,l) 
			_x2 = _x2 * _x2
			d0_2 = CheckREDC(rmod,n,n_,_x2,l) 
			return ((d0_1, d0_2), (d1_1, d1_2))
		else:			
			if i == '0':
				_x2 = _x1 * _x2
				_x2 = REDC(rmod,n,n_,_x2,l) 
				_x1 = _x1 * _x1
				_x1 = REDC(rmod,n,n_,_x1,l) 
			else:
				_x1 = _x1 * _x2
				_x1 = REDC(rmod,n,n_,_x1,l) 
				_x2 = _x2 * _x2
				_x2 = REDC(rmod,n,n_,_x2,l) 
			c -= 1
			
def CheckDivExp_firstbit(mes,e,n,bit):
	r = findR(n)[1] 
	rmod = r - 1  	
	l = findR(n)[0] 
	n_ = - modinv(n,r) & rmod 
	r2 = (r << l) % n 
	d1 = CheckREDC(rmod,n,n_,mes*r2,l) 
	_x1 = REDC(rmod,n,n_,mes*r2,l) 
	_x2 = _x1 * _x1
	d2 = CheckREDC(rmod,n,n_,_x2,l)
	return d1, d2
	
def CalcDiv(bit0, bit1):
	di = lambda p0, p1: abs(p0[0]-p1[0]) + abs(p0[1]-p1[1])
	return list(map(di, bit0, bit1))

def Padding8 (n): 
	hex_n = hex(n).rstrip("L").lstrip("0x")
	# print("%s\n" % hex_n)
	padding = 8 - (len(hex_n) % 8);
	for i in range(padding):
		hex_n = "0" + hex_n;
	return hex_n;
		
def IsNoDiv (num, mod, e, bit, r1, r2):					
		d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
		if CalcDiv(d1,d2)[0] == CalcDiv(d1,d2)[1] == 0:				
			return 1
		else:
			return 0

def GenBranchCombo (num, mod, e, bit, a, b, c, d, f): 	
	for i in range(num):
		while(True):
			r1 = random.randint(2, mod)
			d1 = CheckDivExp(r1, e, mod, bit)
			if d1[0][0] == a and d1[0][1] == b and d1[1][0] == c and d1[1][1] == d:
				# print("%s, %s" % ( Padding8(r1), str(d1) ) )	
				# print("%d\n" % r1)			
				f.write("%s\n" % Padding8(r1))
				break

def FindDiv (num, mod, e, bit): 
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] >= 1 and CalcDiv(d1,d2)[1] == 0:				
				print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break

def FindNoDiv (num, mod, e, bit): 
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] == CalcDiv(d1,d2)[1] == 0:				
				print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break

def FindDiv_ex (num, mod, e, bit): 
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] >= 1 and CalcDiv(d1,d2)[1] == 0:				
				print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break

def FindNoDiv_ex (num, mod, e, bit): 
	for i in range(num):
		while(True):
			r1, r2 = random.randint(2, mod), random.randint(2, mod)
			d1, d2 = CheckDivExp(r1, e, mod, bit), CheckDivExp(r2, e, mod, bit)
			if CalcDiv(d1,d2)[0] == CalcDiv(d1,d2)[1] == 0:				
				print ("%d, %d, %s, %s" % (r1, r2, str(d1), str(d2)))
				break


random.seed(time.time())
p = 32416189867
q = 32416189909
n = p*q 
phi = (p-1)*(q-1)
n_lambda = phi // egcd(p-1, q-1)[0] 
e = 5
d = modinv(e, n_lambda) #67 bits

# ############ first bit
# f= open("1branchcombo0000_1.txt","w+")
# GenBranchCombo( 1001, n, e, 1, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("2branchcombo0000_1.txt","w+")
# GenBranchCombo( 1001, n, e, 1, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("3branchcombo0000_1.txt","w+")
# GenBranchCombo( 1001, n, e, 1, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("4branchcombo0100_1.txt","w+")
# GenBranchCombo( 1001, n, e, 1, 0, 1, 0, 0, f)
# f.close()
# 
# ############ second bit
# f= open("1branchcombo0000_0.txt","w+")
# GenBranchCombo( 1001, n, e, 0, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("2branchcombo0000_0.txt","w+")
# GenBranchCombo( 1001, n, e, 0, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("3branchcombo0000_0.txt","w+")
# GenBranchCombo( 1001, n, e, 0, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("4branchcombo0100_0.txt","w+")
# GenBranchCombo( 1001, n, e, 0, 0, 1, 0, 0, f)
# f.close()



############ first bit
f= open("1branchcombo0000_65.txt","w+")
GenBranchCombo( 1001, n, d, 65, 0, 0, 0, 0, f)
f.close()

f= open("2branchcombo0000_65.txt","w+")
GenBranchCombo( 1001, n, d, 65, 0, 0, 0, 0, f)
f.close()

f= open("3branchcombo0000_65.txt","w+")
GenBranchCombo( 1001, n, d, 65, 0, 0, 0, 0, f)
f.close()

f= open("4branchcombo0100_65.txt","w+")
GenBranchCombo( 1001, n, d, 65, 0, 1, 0, 0, f)
f.close()

############ second bit
f= open("1branchcombo0000_64.txt","w+")
GenBranchCombo( 1001, n, d, 64, 0, 0, 0, 0, f)
f.close()

f= open("2branchcombo0000_64.txt","w+")
GenBranchCombo( 1001, n, d, 64, 0, 0, 0, 0, f)
f.close()

f= open("3branchcombo0000_64.txt","w+")
GenBranchCombo( 1001, n, d, 64, 0, 0, 0, 0, f)
f.close()

f= open("4branchcombo0100_64.txt","w+")
GenBranchCombo( 1001, n, d, 64, 0, 1, 0, 0, f)
f.close()

# PTX, more samples, don't use combo, stop at divergent bit (shorter bits), making preceding bits non divergent, cahces, 32 threads.




