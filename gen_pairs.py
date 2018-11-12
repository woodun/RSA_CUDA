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
			
def CheckDivExp(mes1, mes2, e, n, bit, check_pre):
	r = findR(n)[1] 
	rmod = r - 1  	
	l = findR(n)[0] 
	n_ = - modinv(n,r) & rmod 
	r2 = (r << l) % n 
	
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
	
	e_b = bits(e)
	if bit > len(e_b) - 2 :
		print ("Wrong bit!")
		exit(1)
	c = len(e_b) - 2
	
	for i in e_b[1:]:
		
		if check_pre == 1 and ( s1_1 != s1_2 or s2_1 != s2_2 ): #previous bits are all convergent no matter it is 0 or 1
			return 0
		
		if bit == c:			
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
					#print ("debug0\n")
					return 0
				else: #diverge for bit 0 and converge for bit 1
					#print ("debug4\n")
					return 4
			elif d1_s1_1 != d1_s1_2 or d1_s2_1 != d1_s2_2: #converge for bit 0, diverge for bit 1
				#print ("debug1\n")
				return 1
			else: #converge for bit 0 and converge for bit 1
				#print ("debug\n")
				return 2			
		else:
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
			c -= 1
			
def Padding8 (n): 
	hex_n = hex(n).rstrip("L").lstrip("0x")
	# print("%s\n" % hex_n)
	padding = 8 - (len(hex_n) % 8);
	for i in range(padding):
		hex_n = "0" + hex_n;
	return hex_n;
			
def FindPairs (num, mod, e, bit, f1, f2, f3, check_pre): 
	bit1_div_num = num
	nondiv_num = num
	bit0_div_num = num
	while(True):
		r1, r2 = random.randint(2, mod), random.randint(2, mod)
		div_con = CheckDivExp(r1, r2, e, mod, bit, check_pre)	
		if div_con == 1 and bit1_div_num > 0:				
			f1.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			bit1_div_num-=1
		if div_con == 2 and nondiv_num > 0:				
			f2.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			nondiv_num-=1
		if div_con == 4 and bit0_div_num > 0:				
			f3.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
			bit0_div_num-=1
		if bit1_div_num == 0 and nondiv_num == 0 and bit0_div_num == 0:
			break

random.seed(time.time())
p = 32416189867
q = 32416189909
n = p * q
phi = (p - 1) * (q - 1)
n_lambda = phi // egcd(p-1, q-1)[0] 
e = 5
d = modinv(e, n_lambda) #67 bits


# PTX, more samples, don't use combo, stop at divergent bit (shorter bits), making preceding bits non divergent, cahces, 32 threads. If check pre all convergent not working, try all divergent or 1/3 divergent.
# f1 = open("divpairs_nopre1.txt","w+")
# f2 = open("nondivpairs_nopre1.txt","w+")
# FindPairs (10000, n, e, 1, f1, f2, 0)
# f1.close()
# f2.close()
#  
# f1 = open("divpairs_nopre0.txt","w+")
# f2 = open("nondivpairs_nopre0.txt","w+")
# FindPairs (10000, n, e, 0, f1, f2, 0)
# f1.close()
# f2.close()
# 
# f1 = open("divpairs_pre1.txt","w+")
# f2 = open("nondivpairs_pre1.txt","w+")
# FindPairs (10000, n, e, 1, f1, f2, 1)
# f1.close()
# f2.close()
#  
# f1 = open("divpairs_pre0.txt","w+")
# f2 = open("nondivpairs_pre0.txt","w+")
# FindPairs (10000, n, e, 0, f1, f2, 1)
# f1.close()
# f2.close()


# f1 = open("divpairs_nopre65.txt","w+")
# f2 = open("nondivpairs_nopre65.txt","w+")
# FindPairs (10000, n, d, 65, f1, f2, 0)
# f1.close()
# f2.close()
#  
# f1 = open("divpairs_nopre64.txt","w+")
# f2 = open("nondivpairs_nopre64.txt","w+")
# FindPairs (10000, n, d, 64, f1, f2, 0)
# f1.close()
# f2.close()
# 
# f1 = open("divpairs_pre65.txt","w+")
# f2 = open("nondivpairs_pre65.txt","w+")
# FindPairs (10000, n, d, 65, f1, f2, 1)
# f1.close()
# f2.close()
 
# f1 = open("bit1divpairs_pre64.txt","w+")
# f2 = open("nondivpairs_pre64.txt","w+")
# f3 = open("bit0divpairs_pre64.txt","w+")
# FindPairs (10000, n, d, 64, f1, f2, f3, 1)
# f1.close()
# f2.close()
# f3.close()

f1 = open("bit1divpairs_pre64.txt","w+",1)
f2 = open("nondivpairs_pre64.txt","w+",1)
f3 = open("bit0divpairs_pre64.txt","w+",1)
FindPairs (10, n, d, 31, f1, f2, f3, 0)
f1.close()
f2.close()
f3.close()


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
# f= open("1branchcombo0000_65.txt","w+")
# GenBranchCombo( 1001, n, d, 65, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("2branchcombo0000_65.txt","w+")
# GenBranchCombo( 1001, n, d, 65, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("3branchcombo0000_65.txt","w+")
# GenBranchCombo( 1001, n, d, 65, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("4branchcombo0100_65.txt","w+")
# GenBranchCombo( 1001, n, d, 65, 0, 1, 0, 0, f)
# f.close()
# 
# ############ second bit
# f= open("1branchcombo0000_64.txt","w+")
# GenBranchCombo( 1001, n, d, 64, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("2branchcombo0000_64.txt","w+")
# GenBranchCombo( 1001, n, d, 64, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("3branchcombo0000_64.txt","w+")
# GenBranchCombo( 1001, n, d, 64, 0, 0, 0, 0, f)
# f.close()
# 
# f= open("4branchcombo0100_64.txt","w+")
# GenBranchCombo( 1001, n, d, 64, 0, 1, 0, 0, f)
# f.close()






