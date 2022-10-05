'''
Created on 5 mars 2017

@author: garzol
'''

'''
-----------------------------
Miller-rabin primality test
-----------------------------
'''

from random import randint 
from gmpy2 import gcd

def isWitness(possibleWitness, p, exponent, remainder):
    possibleWitness = pow(possibleWitness, remainder, p)

    if possibleWitness == 1:
        return False
    
    prev = possibleWitness
    for _ in range(exponent):
        possibleWitness = pow(possibleWitness, 2, p)
        if possibleWitness == 1:
            g = gcd(prev+1,p)
            if g==1 or g ==p:
                return False
            return g
        prev = possibleWitness
    return True

def probablyPrime(p, accuracy=20):
    if p == 2 or p == 3: return True
    if p < 2: return False
    
    exponent, remainder = 0, p-1
    while remainder & 1 == 0:
        remainder >>= 1
        exponent += 1

    for _ in range(accuracy):
        possibleWitness = randint(2, p - 2)
        k = isWitness(possibleWitness, p, exponent, remainder)
        if k == True:
            return False
        if k:
            return k
    return True



    