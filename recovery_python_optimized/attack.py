from __future__ import division

import random, time, sys, re, os, subprocess

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
            
def CheckDivExp(mes1, mes2, e, n, n_, rsquare, rmod, l):
    
    s1_1 = CheckREDC(rmod, n, n_, mes1 * rsquare, l)
    s1_2 = CheckREDC(rmod, n, n_, mes2 * rsquare, l)    
    _x1_1 = REDC(rmod, n, n_, mes1 * rsquare, l) 
    _x1_2 = REDC(rmod, n, n_, mes2 * rsquare, l)
    
    _x2_1 = _x1_1 * _x1_1
    _x2_2 = _x1_2 * _x1_2
    s2_1 = CheckREDC(rmod, n, n_, _x2_1, l)
    s2_2 = CheckREDC(rmod, n, n_, _x2_2, l)
    _x2_1 = REDC(rmod, n, n_, _x2_1, l)
    _x2_2 = REDC(rmod, n, n_, _x2_2, l)
    
    e_b = bits(e)   
    
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
    
    # we do not use (1 1)(0 0) because it's never gonna happen
    if (d0_s1_1 != d0_s1_2 and d0_s2_1 == d0_s2_2) or (d0_s1_1 == d0_s1_2 and d0_s2_1 != d0_s2_2): #diverge for bit 0 (1 0) or (0 1)
        if (d1_s1_1 != d1_s1_2 and d1_s2_1 == d1_s2_2) or (d1_s1_1 == d1_s1_2 and d1_s2_1 != d1_s2_2): #diverge for bit 0, diverge for bit 1 (1 0) or (0 1)
#             print ("debug3\n")
            return 3
        elif d1_s1_1 == d1_s1_2 and d1_s2_1 == d1_s2_2: #diverge for bit 0, converge for bit 1 (0 0)
#             print ("debug4\n")
            return 4
        else:
            return 0
    elif d0_s1_1 == d0_s1_2 and d0_s2_1 == d0_s2_2: #converge for bit 0 (0 0)
        if (d1_s1_1 != d1_s1_2 and d1_s2_1 == d1_s2_2) or (d1_s1_1 == d1_s1_2 and d1_s2_1 != d1_s2_2): #converge for bit 0, diverge for bit 1 (1 0) or (0 1)
#             print ("debug1\n")
            return 1
        elif d1_s1_1 == d1_s1_2 and d1_s2_1 == d1_s2_2: #converge for bit 0, converge for bit 1 (0 0)
#             print ("debug2\n")
            return 2
        else:
            return 0
    else:
        return 0   
            
def Padding8 (n): 
    hex_n = hex(n).rstrip("L").lstrip("0x")
    # print("%s\n" % hex_n)
    padding = 8 - (len(hex_n) % 8);
    for i in range(padding):
        hex_n = "0" + hex_n;
    return hex_n;
            
def FindPairs (num, mod, e, n_, rsquare, rmod, l, f1, f4): 
    bit1_div_num = num #0 1
    bit0_div_num = num #1 0
    while(True):
        mes1, mes2 = random.randint(2, mod), random.randint(2, mod)
        div_con = CheckDivExp(mes1, mes2, e, mod, n_, rsquare, rmod, l)    
        if div_con == 1 and bit1_div_num > 0:                
            f1.write("%s\n%s\n" % (Padding8(mes1), Padding8(mes2) ) )
            bit1_div_num-=1
        if div_con == 4 and bit0_div_num > 0:                
            f4.write("%s\n%s\n" % (Padding8(mes1), Padding8(mes2) ) )
            bit0_div_num-=1
        if bit1_div_num == 0 and bit0_div_num == 0: # all               
            return 1
          
time1 = time.time()      
print(time1)    
random.seed(time1)
p = 32416189867
q = 32416189909

p = 29996224275497
q = 29996224275821

n = p * q
r = findR(n)[1] 
rmod = r - 1      
l = findR(n)[0] 
n_ = - modinv(n,r) & rmod 
r2 = (r << l) % n

phi = (p-1)*(q-1)
n_lambda = phi // egcd(p-1, q-1)[0]
e = 7
_d = modinv(e, phi)
d = modinv(e, n_lambda)
print( "CUDA inputs: hex(n):%s, hex(N_):%s, hex(R2):%s" % ( Padding8 (n), Padding8 (n_), Padding8 (r2) ) )
print( "CUDA inputs: bits(e):%s, bits(d):%s, len(bits(d)):%s, RL:%d\n" % ( bits(e), bits(d), len(bits(d)), l ) )

#100 zeros
#10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
exit()
_key = "1011011001001001010011110110010101010111001010110101111000111100001"
#key = "1000100010110110111110111000110000000001011000001000011010101101000101"
key = bits(d)

current_bits = 1
eob = 0
temp = "0"
bit_count = 0

print("current bits: " + bits(current_bits))
while(eob == 0 ):    

    f1 = open("bit1divpairs_pre.txt","w+")
    f4 = open("bit0divpairs_pre.txt","w+")

    FindPairs (int(sys.argv[1]), n, current_bits, n_, r2, rmod, l, f1, f4)
    f1.close()
    f4.close()
    
    #./main bit0divpairs_pre.txt 1000
    sum1 = subprocess.check_output(["./main", "bit1divpairs_pre.txt", sys.argv[1], "bit1divpairs_out.txt"]) # greater means 1
    print(sum1)
    sum4 = subprocess.check_output(["./main", "bit0divpairs_pre.txt", sys.argv[1], "bit0divpairs_out.txt"]) # greater means 0
    print(sum4)

    diff5 = int(sum1) - int(sum4) # greater means 1
    print(diff5)    

    if abs(diff5) < 1000 : # noise filter
        print("bit not accepted. current bits: " + bits(current_bits))
        continue

    if diff5 > 0 : #bit is 1
        print("bit accepted: 1")
        temp = bits(current_bits) + "1"
    else : #bit is 0
        print("bit accepted: 0")
        temp = bits(current_bits) + "0"

    bit_count+=1
    current_bits = int(temp, 2)
    print("current bits: " + bits(current_bits))
       
    if bits(current_bits)[bit_count] != key[bit_count] : # length of key
        print("wrong key!");
        break
    
    if len(bits(current_bits)) == 70 : # length of key
        break


















