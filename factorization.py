# -*-coding:UTF-8 -*

'''
Created on 9 févr. 2017

@author: garzol
'''

from _functools import partial
from math import log, sqrt
from random import randint
from time import time

from gmpy2 import gcd, mpz, gcdext

from curves import EllCurveEdwards
from miller_rabin import probablyPrime
from parallel import parallelize
from polys import prod_multi_eval, poly_from_roots, diff_tab_dickson, precompute, \
    inverse, CachedFFT


def sieve(n):
    """ Returns  a list of primes < n """
    sieve = [True] * n
    for i in range(3,int(n**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((n-i*i-1)//(2*i)+1)
    return [2] + [i for i in range(3,n,2) if sieve[i]]

def step1(curve,Q,n,B1,primes):
    time1 = time()
    for p in primes:
        if p>B1:
            print("step 1 : ",(time()-time1)*1000, "ms")
            return Q
        w = int(log(int(B1),p))
        for _ in range(w):
            Q = curve.mul(p, Q)
        g = gcd(Q[0],n)
        if g != 1 : 
            return g
    print("step 1 : ",(time()-time1)*1000, "ms")
    return Q

def step2(curve, Q, B1, B2, deg_dickson):
    # x1.4 to "simulate" stage 2 blocks with k=2. We are very lazy
    d = 6*int(sqrt(B2)//(6*1.4))
    vec = diff_tab_dickson(int((B1//d)*d), deg_dickson, d)
    time2 = time()
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
    time3 = time()
    print("Computing roots of G took ",(time3-time2)*1000,"ms")
    vec = diff_tab_dickson(1,deg_dickson,6)
    time4 = time()
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
    time5 = time()
    print("Computing roots of F took ",(time5-time4)*1000,"ms")
    polys = precompute(curve.n,taus)
    B =polys.val
    time6 = time()
    print("Building F from its roots took ",(time6-time5)*1000,"ms")
    Br_inv = inverse(B[::-1],len(B)-1,curve.n)
    time7 = time()
    print("Computing 1/F took ",(time7-time6)*1000,"ms")
    polynomial = poly_from_roots(sigmas,B,Br_inv,curve.n)
    time8 = time()
    print("Building G from its roots took ",(time8-time7)*1000,"ms")
    pme = prod_multi_eval(polynomial, curve.n,polys,Br_inv,taus)
    g = gcd(pme,curve.n)
    return g


def try_factor(n,B,primes,B2):
    timethen = time()
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
    time1 = time()
    #print("step 1 took ",(time1-timethen)*1000," ms")
    WCurve = curve.toWeierstrass(Q)
    if not isinstance(WCurve,list):
        return WCurve
    ret = step2(WCurve[0],WCurve[1],B,B2,30)
    #print("step 2 took ",(time()-time1)*1000," ms")
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
            n//=k
            while n%k==0:
                n//=k
            contenders.append(n)
    return champs
    
