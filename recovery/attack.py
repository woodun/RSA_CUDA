import re
import sys
import os
import subprocess

def kernel_Launch(MesList): # pass randomly generated MesList and run RSA executable. Returns time taken.
    returnValue = subprocess.call( " path to the executable ", MesList)
    return returnValue

def get_Pair_Launch(f1, f2): # call get_pair.py in order to generate input pair
    subprocess.call(" path to the get_pair.py ", f1, f2) # need to edit get_pair.py to return file address


encryptionBit = []
discrepancy = 1000 # subject to change

i = 0
dBitLength = 64
f1 = open("divpairs_pre0.txt", "w+")
f2 = open("nondivpairs_pre0.txt", "w+")
while i<dBitLength:


    get_Pair_Launch(f1, f2) # need to check whether it's still compatible with check_pre

    t_uniform = kernel_Launch(f2) # returns time taken for given input set.
    t_diverse = kernel_Launch(f1)

    if (t_diverge - t_uniform) < discrepancy:
        encryptionBit.append(1)
    else:
        encryptionBit.append(0)

    i += 1

f1.close()
f2.close()
print(encryptionBit)
