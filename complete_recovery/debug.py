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
            

def Exp1(mes1, mes2, e, n, n_, r2, rmod, l, check_pre, div_num): # div_num is a relaxed condition, otherwise cannot reach very far bits. Also the non-bothdiv combo can always be reached.
    
    print("mes1: ", hex(mes1))
    print("mes2: ", hex(mes2))
    
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
    
    div_count = 0
     
    print("mes1: ", hex(_x1_1), hex(_x2_1))
    print("mes2: ", hex(_x1_2), hex(_x2_2))
    print(s1_1, s1_2, s2_1, s2_2)
    if s1_1 != s1_2 :
        div_count+=1
    if s2_1 != s2_2 :
        div_count+=1    
#         if div_count != 1 : #same divergence pattern
#             return 0
#         div_count = 0
    
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
            
        print("mes1: ", hex(_x1_1), hex(_x2_1))
        print("mes2: ", hex(_x1_2), hex(_x2_2))
        print(s1_1, s1_2, s2_1, s2_2)
        if s1_1 != s1_2 :
            div_count+=1
        if s2_1 != s2_2 :
            div_count+=1
#             if div_count != 1 : #same divergence pattern
#                 return 0
#             div_count = 0
        
    s1 = CheckREDC(rmod,n,n_,_x1_1,l)
    s2 = CheckREDC(rmod,n,n_,_x1_2,l)
    _x1_1 = REDC(rmod,n,n_,_x1_1,l)    
    _x1_2 = REDC(rmod,n,n_,_x1_2,l)
           
    print("mes1: ", hex(_x1_1))
    print("mes2: ", hex(_x1_2))
    print(s1, s2)
    if s1 != s2 :
        div_count+=1
        
    print(div_count)
            
    return div_count
            
            
def Exp(mes1, mes2, e, n, n_, r2, rmod, l, check_pre, div_num): # div_num is a relaxed condition, otherwise cannot reach very far bits. Also the non-bothdiv combo can always be reached.
    
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
    
    div_count = 0
     

    if s1_1 != s1_2 :
        div_count+=1
    if s2_1 != s2_2 :
        div_count+=1    
#         if div_count != 1 : #same divergence pattern
#             return 0
#         div_count = 0
    
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
            

        if s1_1 != s1_2 :
            div_count+=1
        if s2_1 != s2_2 :
            div_count+=1
#             if div_count != 1 : #same divergence pattern
#                 return 0
#             div_count = 0
            
    s1 = CheckREDC(rmod,n,n_,_x1_1,l)
    s2 = CheckREDC(rmod,n,n_,_x1_2,l)
       
    if s1 != s2 :
        div_count+=1
            
    return div_count
    
    
def CheckDivExp(mes1, mes2, e, n, n_, r2, rmod, l, check_pre, div_num): # div_num is a relaxed condition, otherwise cannot reach very far bits. Also the non-bothdiv combo can always be reached.
    
    print("mes1: ", hex(mes1))
    print("mes2: ", hex(mes2))
    
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
    
    print("mes1: ", hex(_x1_1), hex(_x2_1))
    print("mes2: ", hex(_x1_2), hex(_x2_2))
    
    div_count = 0
     

    if s1_1 != s1_2 :
        div_count+=1
    if s2_1 != s2_2 :
        div_count+=1    
#         if div_count != 1 : #same divergence pattern
#             return 0
#         div_count = 0
    
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
            

        if s1_1 != s1_2 :
            div_count+=1
        if s2_1 != s2_2 :
            div_count+=1
#             if div_count != 1 : #same divergence pattern
#                 return 0
#             div_count = 0
        print("mes1: ", hex(_x1_1), hex(_x2_1))
        print("mes2: ", hex(_x1_2), hex(_x2_2))
            
    if check_pre == 1 and div_count != div_num : #total divergence number
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
            
def FindPairs (num, mod, e, n_, r2, rmod, l, check_pre, div_num, d): 
    
    bit1_div_sum = 0 #0 1
    nondiv_sum = 0 #0 0
    bothdiv_sum = 0 #1 1
    bit0_div_sum = 0 #1 0
    
    bit1_div_num = num #0 1
    nondiv_num = num #0 0
    bothdiv_num = num #1 1
    bit0_div_num = num #1 0
    while(True):
        mes1, mes2 = random.randint(2, mod), random.randint(2, mod)

        div_con = CheckDivExp(mes1, mes2, e, mod, n_, r2, rmod, l, check_pre, div_num)    
        if div_con == 1 and bit1_div_num > 0:
            bit1_div_sum+=Exp(mes1, mes2, d, mod, n_, r2, rmod, l, check_pre, div_num)           
            bit1_div_num-=1
            print(1)
            print(hex(mes1))
            print(hex(mes2))
            #break
        if div_con == 2 and nondiv_num > 0:
            nondiv_sum+=Exp1(mes1, mes2, d, mod, n_, r2, rmod, l, check_pre, div_num)
            nondiv_num-=1
            print(2)
            print(hex(mes1))
            print(hex(mes2))
            break
        if div_con == 3 and bothdiv_num > 0:
            bothdiv_sum+=Exp(mes1, mes2, d, mod, n_, r2, rmod, l, check_pre, div_num)
            bothdiv_num-=1
            print(3)
            print(hex(mes1))
            print(hex(mes2))
        if div_con == 4 and bit0_div_num > 0:
            bit0_div_sum+=Exp(mes1, mes2, d, mod, n_, r2, rmod, l, check_pre, div_num)
            bit0_div_num-=1
            print(4)
            print(hex(mes1))
            print(hex(mes2))
        if bit1_div_num == 0 and nondiv_num == 0 and bit0_div_num == 0 and bothdiv_num == 0: # all               
            sum1 = bit1_div_sum / num
            sum2 = nondiv_sum / num
            sum3 = bothdiv_sum / num
            sum4 = bit0_div_sum / num
            
            print("#########################################CPU output###########################################")
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
            
            break
          
time1 = time.time()
print(time1)    
# random.seed(time1)
random.seed(0)
p = 32416189867
q = 32416189909
n = p * q
r = findR(n)[1] 
rmod = r - 1      
l = findR(n)[0] 
n_ = - modinv(n,r) & rmod 
r2 = (r << l) % n


phi = (p - 1) * (q - 1)
n_lambda = phi // egcd(p-1, q-1)[0] 
e = 5
_d = modinv(e, n_lambda) #67 bits
d = modinv(e, phi) #70 bits


current_bits = 2
eob = 0
temp = "0"
bit_count = 0

# vote_0 = 0 # multiple runs to vote?
# vote_1 = 0
key = "1011011001001001010011110110010101010111001010110101111000111100001"
#key = "1000100010110110111110111000110000000001011000001000011010101101000101"
# int(sys.argv[1])
for i in range(1):    
    FindPairs (1, n, current_bits, n_, r2, rmod, l, 0, len(bits(current_bits) ), _d)














