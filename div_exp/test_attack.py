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
    
    div_count = 0
     
    if check_pre == 1: 
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
            
        if check_pre == 1: 
            if s1_1 != s1_2 :
                div_count+=1
            if s2_1 != s2_2 :
                div_count+=1
#             if div_count != 1 : #same divergence pattern
#                 return 0
#             div_count = 0
            
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
            
def FindPairs (num, mod, e, n_, r2, rmod, l, f1, f2, f3, f4, check_pre, div_num): 
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
        if div_con == 3 and bothdiv_num > 0:                
            f3.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
            bothdiv_num-=1
        if div_con == 4 and bit0_div_num > 0:                
            f4.write("%s\n%s\n" % (Padding8(r1), Padding8(r2) ) )
            bit0_div_num-=1
#         if bit1_div_num == 0 and bothdiv_num == 0 and bit0_div_num == 0: # no 0 0            
#             return 0    
#         if bothdiv_num == 0 == 0 and nondiv_num == 0 and bit0_div_num == 0: # no 0 1        
#             return 1    
#         if bit1_div_num == 0 and nondiv_num == 0 and bothdiv_num == 0: # no 1 0            
#             return 2       
#         if bit1_div_num == 0 and nondiv_num == 0 and bit0_div_num == 0: # no 1 1                
#             return 3     
        if bit1_div_num == 0 and nondiv_num == 0 and bit0_div_num == 0 and bothdiv_num == 0: # all               
            return 4
          
time1 = time.time()      
print(time1)    
random.seed(time1)
p = 32416189867
q = 32416189909
n = p * q
r = findR(n)[1] 
rmod = r - 1      
l = findR(n)[0] 
n_ = - modinv(n,r) & rmod 
r2 = (r << l) % n

current_bits = 1
eob = 0
temp = "0"
bit_count = 0

# vote_0 = 0 # multiple runs to vote?
# vote_1 = 0
key = "1011011001001001010011110110010101010111001010110101111000111100001"
# int(sys.argv[1])


# try large samples, ad hoc approach. The noise mainly comes from other bits (affected by the pairs generated), not the hardware. Can we do something after or even before the GPU run to know their quality?
print("1")
while(eob == 0 ):    

    f1 = open("bit1divpairs_pre.txt","w+")
    f2 = open("nondivpairs_pre.txt","w+")
    f3 = open("divpairs_pre.txt","w+")
    f4 = open("bit0divpairs_pre.txt","w+")
#    random.seed(1542226726.06)#get a good seed?
#    random.seed(0)#get a good seed?
    FindPairs (20256, n, current_bits, n_, r2, rmod, l, f1, f2, f3, f4, int(sys.argv[1]), len(bits(current_bits) ) )
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    
    #./main bit0divpairs_pre.txt 1000
    out1 = subprocess.check_output(["./main", "bit1divpairs_pre.txt", "20256", "bit1divpairs_out.txt"]) # greater means 1
    sum1 = out1.splitlines()[0]
    print(sum1)
    div1 = out1.splitlines()[1]
    print(div1)        
    print(out1.splitlines()[2])
    break


    out2 = subprocess.check_output(["./main", "nondivpairs_pre.txt", "20256", "nondivpairs_out.txt"])
    sum2 = out2.splitlines()[0]
    print(sum2)
    div2 = out2.splitlines()[1]
    print(div2)
    
    out3 = subprocess.check_output(["./main", "divpairs_pre.txt", "20256", "divpairs_out.txt"])
    sum3 = out3.splitlines()[0]
    print(sum3)
    div3 = out3.splitlines()[1]
    print(div3)
    
    out4 = subprocess.check_output(["./main", "bit0divpairs_pre.txt", "20256", "bit0divpairs_out.txt"]) # greater means 0
    sum4 = out4.splitlines()[0]
    print(sum4)
    div4 = out4.splitlines()[1]
    print(div4)

    diff1 = int(sum1) - int(sum2) # close to zero means 0, greater than zero means 1
    diff2 = int(sum4) - int(sum2) # close to zero means 1, greater than zero means 0
    mean1 = (diff1 + diff2) / 2
    print(diff1,diff2,mean1)
    
    diff3 = int(sum1) - int(sum3) # close to zero means 1, smaller than zero means 0
    diff4 = int(sum4) - int(sum3) # close to zero means 0, smaller than zero means 1
    mean2 = (diff3 + diff4) / 2
    print(diff3,diff4, mean2)
    
    val1 = mean1 + mean2
    sign1 = mean1 * mean2
    print(val1, sign1)
    
    diff5 = int(sum1) - int(sum4) # greater means 1
    diff6 = int(sum3) - int(sum2) # must be greater
    print(diff5,diff6)

#     if ( ( diff3 > -1000 and diff4 > -1000 ) or ( abs(diff3) > 1000 and abs(diff4) > 1000 ) ) :
#         print("bit not accepted.")
#         continue


#     if abs(diff5) < 1000 :
#         print("bit not accepted.")
#         continue

#     if diff6 < 1000 :
#         print("bit not accepted.")
#         continue

#     if diff6 < 1000 :
#         print("bit not accepted.")
#         continue

#     if ( ( diff1 < 1000 and diff2 < 1000 ) or ( abs(diff1) > 1000 and abs(diff2) > 1000 ) ) and ( ( diff3 > -1000 and diff4 > -1000 ) or ( abs(diff3) > 1000 and abs(diff4) > 1000 ) ) :
#         print("bit not accepted.")
#         continue

    if int(sum1) > int(sum4) : #bit is 1
        print("1")
        temp = bits(current_bits) + "1"
    else : #bit is 0
        print("0")
        temp = bits(current_bits) + "0"    

    bit_count+=1
    current_bits = int(temp, 2)
    print("bits: " + bits(current_bits))
       
    if bits(current_bits)[bit_count] != key[bit_count] : # length of key
        print("wrong key!");
        break
    
    if len(bits(current_bits)) == 67 : # length of key
        break


















