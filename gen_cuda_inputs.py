from __future__ import division

import random
import math
import sys

sys.setrecursionlimit(1000000)

def rabinMiller(n):
     s = n-1
     t = 0
     while s&1 == 0:
         s = s//2
         t +=1
     k = 0
     while k<128:
         a = random.randrange(2,n-1)
         #a^s is computationally infeasible.  we need a more intelligent approach
         #v = (a**s)%n
         #python's core math module can do modular exponentiation
         v = pow(a,s,n) #where values are (num,exp,mod)
         if v != 1:
             i=0
             while v != (n-1):
                 if i == t-1:
                     return False
                 else:
                     i = i+1
                     v = (v**2)%n
         k+=2
     return True

def isPrime(n):
     #lowPrimes is all primes (sans 2, which is covered by the bitwise and operator)
     #under 1000. taking n modulo each lowPrime allows us to remove a huge chunk
     #of composite numbers from our potential pool without resorting to Rabin-Miller
     lowPrimes =   [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
                   ,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179
                   ,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269
                   ,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367
                   ,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461
                   ,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571
                   ,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661
                   ,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773
                   ,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883
                   ,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997]
     if (n >= 3):
         if (n&1 != 0):
             for p in lowPrimes:
                 if (n == p):
                    return True
                 if (n % p == 0):
                     return False
             return rabinMiller(n)
     return False

def generateLargePrime(k):
     #k is the desired bit length
     r=100*(math.log(k,2)+1) #number of attempts max
     r_ = r
     while r>0:
        #randrange is mersenne twister and is completely deterministic
        #unusable for serious crypto purposes
         n = random.randrange(2**(k-1),2**(k))
         r-=1
         if isPrime(n) == True:
             return n
     return ("Failure after ", r_, " tries.")

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

def Padding8 (n): 
    hex_n = hex(n).rstrip("L").lstrip("0x")
    # print("%s\n" % hex_n)
    padding = 8 - (len(hex_n) % 8);
    for i in range(padding):
        hex_n = "0" + hex_n;
    return hex_n;


#print( generateLargePrime(1025) )
#p = 231316639582577726089468699223295937290637330268991096759452561191101292162773066141786105564705747159539006801962878343715140187565228438889102198142092209449756532130151626648732655217512366948998506268478869375183672716728155662269808823876705338020319444638588631754295440399401078143964374978373077825387
#q = 208846958141762590157353458492266169455379335074657470476383548260190383198784345734596527533880749852805754278083442767349967987521436175699791397967607078449366149933673731416850907861158380810360655942350418233822286661762728105307362692046980122522432587052138497461328896093288687101759473686425239943461

print( generateLargePrime(1025) )

p = 7
q = 5

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
print( "CUDA inputs:" )
print( "hex(N):%s" % ( Padding8 (n) ) )
print( "hex(N_):%s" % ( Padding8 (n_) ) )
print( "hex(R2):%s" % ( Padding8 (r2) ) )
print( "bits(e):%s" % ( bits(e) ) )
print( "bits(d):%s" % ( bits(d) ) )
print( "len(bits(d)):%s" % ( len(bits(d)) ) )
print( "RL:%d\n" % ( l ) )

#qsub -I -l nodes=1:hima:p100:ppn=1 -l walltime=1:00:00













