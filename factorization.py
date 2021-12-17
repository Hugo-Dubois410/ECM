# -*-coding:UTF-8 -*

'''
Created on 9 févr. 2017

@author: garzol
'''

from _functools import partial
from math import log, sqrt
from random import randint
from time import time

from curves import EllCurveEdwards
from gmpy2 import gcd, mpz, gcdext
from miller_rabin import probablyPrime
from parallel import parallelize
from polys import prod_multi_eval, poly_from_roots, diff_tab_dickson


def sieve(n):
    """ Returns  a list of primes < n """
    sieve = [True] * n
    for i in range(3,int(n**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((n-i*i-1)//(2*i)+1)
    return [2] + [i for i in range(3,n,2) if sieve[i]]


def step1(curve,Q,n,B1,primes):
    startt = time()
    for p in primes:
        if p>B1:
            print("step1 took",1000*(time()-startt),"ms")
            return Q
        w = int(log(int(B1),p))
        for _ in range(w):
            Q = curve.mul(p, Q)
        g = gcd(Q[0],n)
        if g != 1 : 
            return g
    print("step1 took",1000*(time()-startt),"ms")
    return Q

def step2(curve, Q, B1, B2, deg_dickson):
    tstart=time()
    d = 6*int(sqrt(B2)//6)
    vec = diff_tab_dickson(int((B1//d)*d), deg_dickson, d)
    for i in range(len(vec)):
        t = curve.mul(vec[i],Q)
        if t[2]>1:
            return t[2]
        vec[i] = [t[0],t[1]]
    sigmas = [vec[0][0]]
    for _ in range(int((B2-B1)//d+1)):
        curve.add_parallel(vec)
        if not isinstance(vec, list):
            return vec
        sigmas.append(vec[0][0])
    vec = diff_tab_dickson(1,deg_dickson,6)
    for i in range(len(vec)):
        t = curve.mul(vec[i],Q)
        if t[2]>1:
            return t[2]
        vec[i] = [t[0],t[1]]
    taus = [vec[0][0]]
    for j in range(7,d,6):
        curve.add_parallel(vec)
        if not isinstance(vec, list):
            return vec
        if gcd(j,d) == 1: 
            taus.append(vec[0][0])
    polynomial = poly_from_roots(sigmas,curve.n)
    g = gcd(prod_multi_eval(polynomial, curve.n, taus),curve.n)
    print("step2 took",1000*(time()-tstart),"ms")
    return g



def try_factor(n,B,primes,B2):
    x ,y = mpz(randint(0,n-1)), mpz(randint(0,n-1))
    g1, a,_ = gcdext(x,n)
    g2,b,_  = gcdext(y,n)
    if(g1 !=1):
        return g1
    if(g2 !=1):
        return g2
    d = (a*a-b*b*(a*a+1))%n
    P = [x,y,1,(x*y)%n]
    curve = EllCurveEdwards(d,n)
    Q = step1(curve,P,n,B,primes)
    if not isinstance(Q,list):
        return Q
    WCurve = curve.toWeierstrass(Q)
    if not isinstance(WCurve,list):
        return WCurve
    ret = step2(WCurve[0],WCurve[1],B,B2,30)
    return ret


def factor(n,primes,nProc,B1,B2):
    champs = set()
    for p in primes:
        if p>1000:
            break
        if n%p ==0:
            n //= p
            champs.add(p)
        while n%p == 0:
            n //= p
            
    contenders = [n]
    while contenders:
        n = contenders.pop()
        if n <=1 : continue  
        t = probablyPrime(n)
        if t == True: 
            champs.add(int(n))
        elif t : 
            contenders.append(t)
            contenders.append(n//t)
        else:
            k = 1
            while k<2 or k==n:
                print("Trying",nProc,"curves")
                if nProc>1 : k = parallelize(partial(try_factor,n,B1,primes,B2), nProc,n)
                else :  k = try_factor(n,B1,primes,B2)
            print("found a factor:",k,"\n")
            contenders.append(k)
            contenders.append(n//k)
    return champs
    
