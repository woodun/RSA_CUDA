from __future__ import division

import random
import math
import sys

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


print( generateLargePrime(513) )

p = 24459629785364283833111470410697016767033653211208122876885582413633536633041146086554558328286062725490645509691783941884323722081206394402126328961429611
q = 23919455913170500015755783799226625762223555612017221740308873996331017696448333058981571497700037939790397009985964455331321308700955400852467790892848337

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















