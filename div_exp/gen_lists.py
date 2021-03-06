#!/usr/bin/python

import random, time, sys

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
			
def CheckDivExp(mes1, mes2, e, n, n_, r2, rmod, l, check_pre, div_num): # div_num is a relaxed condition, otherwise cannot reach very far bits. Also the non-bothdiv combo can always be reached.
	
	s1_1 = CheckREDC(rmod, n, n_, mes1 * r2, l)
	s1_2 = CheckREDC(rmod, n, n_, mes2 * r2, l)	
	_x1_1 = REDC(rmod, n, n_, mes1 * r2, l) 
	_x1_2 = REDC(rmod, n, n_, mes2 * r2, l)
	
	_x2_1 = _x1_1 * _x1_1
	_x2_2 = _x1_2 * _x1_2
	s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
	s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	_x2_1 = REDC(rmod, n, n_, _x2_1, l)
	_x2_2 = REDC(rmod, n, n_, _x2_2, l)
	
	if check_pre == 1: 
		if s1_1 != s1_2 :
			div_count+=1
		if s2_1 != s2_2 :
			div_count+=1
	
	e_b = bits(e)	
	
	div_count = 0
	
	for i in e_b[1:]:		
		if i == '0':
			_x2_1 = _x1_1 * _x2_1
			s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			_x2_1 = REDC(rmod, n, n_, _x2_1, l) 
			_x1_1 = _x1_1 * _x1_1
			s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			_x1_1 = REDC(rmod, n, n_, _x1_1, l) 
			
			_x2_2 = _x1_2 * _x2_2
			s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			_x2_2 = REDC(rmod, n, n_, _x2_2, l) 
			_x1_2 = _x1_2 * _x1_2
			s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			_x1_2 = REDC(rmod, n, n_, _x1_2, l) 				
		else:
			_x1_1 = _x1_1 * _x2_1
			s1_1 = CheckREDC(rmod, n, n_, _x1_1, l)
			_x1_1 = REDC(rmod, n, n_, _x1_1, l) 
			_x2_1 = _x2_1 * _x2_1
			s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
			_x2_1 = REDC(rmod, n, n_, _x2_1, l)
			
			_x1_2 = _x1_2 * _x2_2
			s1_2 = CheckREDC(rmod, n, n_, _x1_2, l)
			_x1_2 = REDC(rmod, n, n_, _x1_2, l) 
			_x2_2 = _x2_2 * _x2_2
			s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
			_x2_2 = REDC(rmod, n, n_, _x2_2, l)		
			
		if check_pre == 1: 
			if s1_1 != s1_2 :
				div_count+=1
			if s2_1 != s2_2 :
				div_count+=1
			
	if check_pre == 1 and div_count != div_num :
		return 0
	
	_x1_1_temp = _x1_1
	_x2_1_temp = _x2_1
	_x1_2_temp = _x1_2
	_x2_2_temp = _x2_2

	#simulate exp bit 0
	_x2_1 = _x1_1 * _x2_1
	d0_s1_1 = CheckREDC(rmod, n, n_, _x2_1, l) 
	_x1_1 = _x1_1 * _x1_1
	d0_s2_1 = CheckREDC(rmod, n, n_, _x1_1 ,l)
	
	_x2_2 = _x1_2 * _x2_2
	d0_s1_2 = CheckREDC(rmod, n, n_, _x2_2, l) 
	_x1_2 = _x1_2 * _x1_2
	d0_s2_2 = CheckREDC(rmod, n, n_, _x1_2 ,l) 

	#simulate exp bit 1
	_x1_1 = _x1_1_temp
	_x2_1 = _x2_1_temp
	_x1_2 = _x1_2_temp
	_x2_2 = _x2_2_temp
	
	_x1_1 = _x1_1 * _x2_1
	d1_s1_1 = CheckREDC(rmod, n, n_, _x1_1, l) 
	_x2_1 = _x2_1 * _x2_1
	d1_s2_1 = CheckREDC(rmod, n, n_, _x2_1, l) 
	
	_x1_2 = _x1_2 * _x2_2
	d1_s1_2 = CheckREDC(rmod, n, n_, _x1_2, l) 
	_x2_2 = _x2_2 * _x2_2
	d1_s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
	
	if d0_s1_1 != d0_s1_2 or d0_s2_1 != d0_s2_2: #diverge for bit 0
		if d1_s1_1 != d1_s1_2 or d1_s2_1 != d1_s2_2: #diverge for bit 0 and diverge for bit 1
# 					print ("debug3\n")
			return 3
		else: #diverge for bit 0 and converge for bit 1
# 					print ("debug4\n")
			return 4
	elif d1_s1_1 != d1_s1_2 or d1_s2_1 != d1_s2_2: #converge for bit 0, diverge for bit 1
# 				print ("debug1\n")
		return 1
	else: #converge for bit 0 and converge for bit 1
# 				print ("debug2\n")
		return 2	
			
def Padding8 (n): 
	hex_n = hex(n).rstrip("L").lstrip("0x")
	# print("%s\n" % hex_n)
	padding = 8 - (len(hex_n) % 8);
	for i in range(padding):
		hex_n = "0" + hex_n;
	return hex_n;
			
def FindPairs (num, mod, e, n_, r2, rmod, l, f1, f2, f4, check_pre, div_num): 
	bit1_div_num = num #0 1
	nondiv_num = num #0 0
	bothdiv_num = num #1 1
	bit0_div_num = num #1 0
	while(True):
		r1, r2 = random.randint(2, mod), random.randint(2, mod)
		div_con = CheckDivExp(r1, r2, e, mod, n_, r2, rmod, l, check_pre, div_num)	
		if div_con == 1 and bit1_div_num > 0:				
			f1.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			bit1_div_num-=1
		if div_con == 2 and nondiv_num > 0:				
			f2.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			nondiv_num-=1
		if div_con == 4 and bit0_div_num > 0:				
			f4.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			bit0_div_num-=1			
		if bit1_div_num == 0 and nondiv_num == 0 and bit0_div_num == 0: # no 1 1				
			return 3
					
mes1, mes2, e, n, n_, r2, rmod, l, check_pre, div_num
random.seed(time.time())
p = 32416189867
q = 32416189909
n = p * q
r = findR(n)[1] 
rmod = r - 1  	
l = findR(n)[0] 
n_ = - modinv(n,r) & rmod 
r2 = (r << l) % n


f1 = open("bit1divpairs_pre.txt","w+",1)
f2 = open("nondivpairs_pre.txt","w+",1)
f4 = open("bit0divpairs_pre.txt","w+",1)
FindPairs (sys.argv[1], n, current_bits, n_, r2, rmod, l, f1, f2, f4, 1, len(bits(current_bits) ) )
f1.close()
f2.close()
f4.close()





