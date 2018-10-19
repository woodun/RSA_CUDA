#!/usr/bin/python



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
		#print "_a %d _b %d" % (_a, _b)
	return _a

def toMont(i, R, N):
	return (i * R) % N

def findR(i):
	i_b = "{0:b}".format(i)
	return 2**len(i_b)

def REDC(R,N,N_,T):
	m = ((T % R) * N_) % R
	t = (T + m*N) / R
	if t >= N:
		return t - N
	else:
		return t
	
# def MontExp(mes,e,n):
# 	r = findR(n)
# 	n_ = - modinv(n,r) % r
# 	_a = toMont(1,r,n)
# 	_b = toMont(mes,r,n)
# 	for i in bits(e)[::-1]:
# 		if i == '1':
# 			_a = _a * _b
# 			_a =  REDC(r,n,n_,_a)
# 		_b = _b * _b
# 		_b = REDC(r,n,n_,_b)
# 		#print "_a %d _b %d" % (_a, _b)
# 	#final reduction to convert from Montgomery into normal form
# 	_a = REDC(r,n,n_,_a)
# 	return _a


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

# e_b = bits(e)
# d_b = bits(d)
# n_b = bits(n)

#print "p: %d q: %d n: %d(%s) phi: %d e: %d(%s) d: %d(%s)" % (p,q,n,n_b,phi,e,e_b,d,d_b)


# p = 61
# q = 53

p = 32416189867
q = 32416189909

p = 61
q = 53

n = p*q # = 3233
phi = (p-1)*(q-1) # = 3120
e = 19
# d = 5795

d = modinv(e, phi)



#encrypt:
mes = 1234
c = pow(mes, e, n)
#decrypt:
m1 = pow(c,d,n)
#print "m: %d c: %d m_dec: %d" % (mes,c,m1)
print (mes,c,m1)
 
#square-and-multiply
c = sqm(mes,e,n)
m2 = sqm(c,d,n)
#print "m: %d c: %d m_dec: %d" % (mes,c,m1)
print (mes,c,m2)
# 
# Testing Montgomery multiplication:
#  R = findR(n)
#  a = 137
#  b = 262
#  c = a * b % n
#  am = toMont(a,R,n)
#  bm = toMont(b,R,n)
#  cm = toMont(c,R,n)
#  N_ = - modinv(n,R) % R
#  tmp1 = am*bm
#  tmp2 = REDC(R,n,N_,tmp1)
#  tmp3 = REDC(R,n,N_,tmp2)
# print "R: %d (%s) N_: %d a: %d am: %d b: %d bm: %d c: %d cm:%d" % (R, bits(R), N_, a, am, b, bm, c, cm)
# print "am*bm=%d, reduced: %d, double-reduced: %d " % (tmp1, tmp2, tmp3)
# 
# Exponentation with Montgomery multiplication:
c = MontExp(mes,e,n)
m2 = MontExp(c,d,n)
# print "m: %d c: %d m_dec: %d" % (mes,c,m1)
print (mes,c,m2)


#print(bits(14)[::-1])
# print(findR(14))
#           
# g, x, y = egcd(17, 780)
#           
# print(g, x, y)
#           
# print(x % 780)
#           
# print(egcd(17, 780))
#           
# print(egcd(106, 82))#53 * 41 = 2173
#           
# print(modinv(3,2173))
#           
#n_ = - modinv(17,780) % 780
#           
#print(n_ * 17 % 780)
#           
# print((-17 * n_)%780)
#     
    












