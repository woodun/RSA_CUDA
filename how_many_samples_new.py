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
			
def CheckDivExp(e, n, n_, r2, rmod, l, bit, check_pre, div_num, num): #now gen pair and check pre at the same time
	
	bit1_div_sum = 0 #0 1
	nondiv_sum = 0 #0 0
	bothdiv_sum = 0 #1 1
	bit0_div_sum = 0 #1 0
	
	bit1_div_num = num #0 1
	nondiv_num = num #0 0
	bothdiv_num = num #1 1
	bit0_div_num = num #1 0
	
	e_b = bits(e)
	bit_length = len(e_b) - 2
	if bit > bit_length :
		print ("Wrong bit!")
		exit(1)
	
	while bit1_div_num > 0 or nondiv_num > 0 or bothdiv_num > 0 or bit0_div_num > 0 :
	
		mark = 0
		div_count = 0
		c = bit_length
		
		mes1, mes2 = random.randint(2, n), random.randint(2, n)
		
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
		
		for i in e_b[1:]:
			
			if s1_1 != s1_2 :
				div_count+=1				
			if s2_1 != s2_2 :
				div_count+=1				
# 			if div_count != 1 : #same divergence pattern
# 				return 0
# 			div_count = 0
			
			if bit == c:			
				
				if check_pre == 1 and div_count != div_num : #total divergence number
					break			
				
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
				
				if (d0_s1_1 != d0_s1_2 and d0_s2_1 == d0_s2_2) or (d0_s1_1 == d0_s1_2 and d0_s2_1 != d0_s2_2): #diverge for bit 0 (1 0) or (0 1)
					if (d1_s1_1 != d1_s1_2 and d1_s2_1 == d1_s2_2) or (d1_s1_1 == d1_s1_2 and d1_s2_1 != d1_s2_2): #diverge for bit 0, diverge for bit 1 (1 0) or (0 1)
						#print ("debug3\n")
						if bothdiv_num > 0 :
							bothdiv_num-=1
							mark = 3
						else:
							break
					elif d1_s1_1 == d1_s1_2 and d1_s2_1 == d1_s2_2: #diverge for bit 0, converge for bit 1 (0 0)
						#print ("debug4\n")
						if bit0_div_num > 0 :
							bit0_div_num-=1
							mark = 4
						else:
							break
					else:
						break
				elif d0_s1_1 == d0_s1_2 and d0_s2_1 == d0_s2_2: #converge for bit 0 (0 0)
					if (d1_s1_1 != d1_s1_2 and d1_s2_1 == d1_s2_2) or (d1_s1_1 == d1_s1_2 and d1_s2_1 != d1_s2_2): #converge for bit 0, diverge for bit 1 (1 0) or (0 1)
						#print ("debug1\n")
						if bit1_div_num > 0 :
							bit1_div_num-=1
							mark = 1
						else:
							break
					elif d1_s1_1 == d1_s1_2 and d1_s2_1 == d1_s2_2: #converge for bit 0, converge for bit 1 (0 0)
						#print ("debug2\n")
						if nondiv_num > 0 :
							nondiv_num-=1
							mark = 2
						else:
							break
					else:
						break
				else:
					break 	
			
				#still continue to finish all bits
				_x1_1 = _x1_1_temp
				_x2_1 = _x2_1_temp
				_x1_2 = _x1_2_temp
				_x2_2 = _x2_2_temp
				
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
				
		if mark == 0 :
			continue
			
		s1 = CheckREDC(rmod,n,n_,_x1_1,l)
		s2 = CheckREDC(rmod,n,n_,_x1_2,l)
		
		
		if s1_1 != s1_2 :
			div_count+=1				
		if s2_1 != s2_2 :
			div_count+=1
		if s1 != s2 :
			div_count+=1
				
		if mark == 1 :
			bit1_div_sum+=div_count
		elif mark == 2 :
			nondiv_sum+=div_count
		elif mark == 3 :
			bothdiv_sum+=div_count
		else: # mark == 4
			bit0_div_sum+=div_count			
			
	sum1 = bit1_div_sum / num
	sum2 = nondiv_sum / num
	sum3 = bothdiv_sum / num
	sum4 = bit0_div_sum / num
	print(sum1)
	print(sum2)
	print(sum3)
	print(sum4)
	
	diff1 = sum1 - sum2 # close to zero means 0, greater than zero means 1
	diff2 = sum4 - sum2 # close to zero means 1, greater than zero means 0
	mean1 = (diff1 + diff2) / 2
	print(diff1,diff2,mean1)
    
	diff3 = sum1 - sum3 # close to zero means 1, smaller than zero means 0
	diff4 = sum4 - sum3 # close to zero means 0, smaller than zero means 1
	mean2 = (diff3 + diff4) / 2
	print(diff3,diff4, mean2)
    
	val1 = mean1 + mean2
	sign1 = mean1 * mean2
	print(val1, sign1)
    
	diff5 = sum1 - sum4 # greater means 1
	diff6 = sum3 - sum2 # must be greater
	print(diff5,diff6)


random.seed(time.time())
p = 32416189867
q = 32416189909
n = p * q
phi = (p - 1) * (q - 1)
n_lambda = phi // egcd(p-1, q-1)[0] 
e = 5
_d = modinv(e, n_lambda) #67 bits
d = modinv(e, phi)

r = findR(n)[1] 
rmod = r - 1  	
l = findR(n)[0] 
n_ = - modinv(n,r) & rmod 
r2 = (r << l) % n
	
bit = 68

#key = "1011011001001001010011110110010101010111001010110101111000111100001"
key = "1000100010110110111110111000110000000001011000001000011010101101000101"
print(bits(d))
print("\n") 

start = time.time()

print(key[(69 - bit)])
for i in range(10):
	CheckDivExp(d, n, n_, r2, rmod, l, bit, 0, (69 - bit), 2000)

end = time.time()

print("\n") 

print(end - start) 





















