'''
Created on 5 juin 2017

@author: garzol
'''

from math import log
from time import time

from gmpy2 import mpz, bit_mask, bit_length, powmod


class Node():
    def __init__(self,val,right=None,left=None):
        self.val = val
        self.right= right
        self.left = left

#Make sure the coefficients are < n
def horner_2power(x, polynomial):
    result = mpz(0)
    for i,coefficient in enumerate(polynomial):
        result ^= (coefficient<<(x*i))
    return result


def mul_kronecker(p,q,n):
    d = max(len(p),len(q))
    numbits = bit_length(n*n*d)
    A = horner_2power(numbits,p)
    B = horner_2power(numbits,q)
    R = A*B
    i = len(p)+len(q)-2
    res = [0]*(i+1)
    mask = bit_mask(numbits)
    i = 0
    while R > 0:
        res[i] = (R & mask )% n 
        i += 1
        R >>= numbits
    return res

def minus(p,q,n):
    if len(p)<len(q) : return [(coeff1-coeff2)%n for (coeff1,coeff2) in zip(p+[mpz(0)]*(len(q)-len(p)),q)]
    else:              return [(coeff1-coeff2)%n for (coeff1,coeff2) in zip(p,q+[mpz(0)]*(len(p)-len(q)))]


def inverse(B,t,n):
    if t<=1:
        return [1]
    if t&1: t0 = (t>>1)+1
    else:   t0 = t>>1
    B0 = inverse(B,t0,n)
    m = mul(B0,B[:t],n)[:t]
    m[0]-=1
    m1 = mul(m,B0,n)[:t]
    return minus(B0,m1,n)

def newton_inv(B,t,n):
    r = int(log(t,2))+1
    res = [1]
    for i in range(1,r+1):
        g2 = mul(res,res,n)
        res = minus([(2*k%n) for k in res], mul(B[:(1<<i)], g2, n), n)[:(1<<i)]
    return res[:t]
'''
def newton_schonhage(g,f,s,n):
    if s%2==0:
        m = t = s//2
        eta = 2
    else:
        m = s//2
        t = m+1
        eta = 1
    s1, m1 , t1 = 1<<s, 1<<m , 1<<t
    P_ = []
    Q_ = []
    j = 0
    for i in range(0,s1,m1):
        P_.append(mulX(g[i:i+m1],j,2*m1))
        Q_.append(mulX(f[i:i+m1],j,2*m1))
        j+=eta
    a = DFT(g,t1,2*eta,2*m1,n)
    return []
'''

def rem(p,B,n):
    Br_inv = inverse(B[::-1],len(p)-len(B)+1,n)
    Dr = mul(p[::-1],Br_inv,n)[:len(p)-len(B)+1]
    return minus(p,mul(Dr[::-1],B,n),n)[:len(B)-1]

def rem_from_inv(p,B,Br_inv,n):
    Dr = mul(p[::-1], Br_inv, n)[:len(p)-len(B)+1]
    return minus(p,mul(Dr[::-1],B,n),n)[:len(B)-1]

def precompute(n,val):
    d = len(val)
    if d==1:
        return Node([n-val[0],1])
    elif d==0:
        return None
    else:
        d>>=1
        right = precompute(n,val[d:])
        left = precompute(n,val[:d])
        return Node(mul(right.val,left.val,n),right,left)

        
def prod_multi_eval(p,n,polys,F_inv,val):
    time2 = time()
    F = polys.val
    #inv_F = inverse(F[::-1],len)
    H = rem_from_inv(p, F, F_inv,n)
    time3 = time()
    print("Reducing G mod F took ",(time3-time2)*1000,"ms")
    res = prod_multi_eval_(H, n, val,polys)
    time4 = time()
    print("Computing polyeval took ",(time4-time3)*1000,"ms")
    return res

def prod_multi_eval_(p,n,val,polys):
    d = len(val)
    if d==0:
        return mpz(1)
    if d==1: 
        return p[0]
    d >>= 1
    f0 = rem(p,polys.left.val,n)
    h1 = prod_multi_eval_(f0,n,val[:d],polys.left)
    polys.left = None
    f1 = rem(p,polys.right.val,n)
    h2 = prod_multi_eval_(f1, n, val[d:],polys.right)
    polys.right = None
    return (h1*h2)%n
    
def poly_from_roots(roots,B,Br_inv, n):
    if len(roots) == 1 : return [n-roots[0],1]
    d = len(roots)>>1
    P = poly_from_roots(roots[d:],B,Br_inv ,n)
    if len(P)>=len(B):
        P = rem_from_inv(P, B, Br_inv,n)
    Q = poly_from_roots(roots[:d],B,Br_inv,n)
    if len(Q)>=len(B):
        Q = rem_from_inv(Q,B,Br_inv,n)
    return mul(P,Q,n)

def dickson(n,x):
    prev, cur = x, x*x+2
    if n<=1: return prev
    for _ in range(2,n):
        prev,cur = cur,x*cur+ prev
    return cur

def diff_tab_dickson(begin,n,h):
    ls = list()
    for _ in range(n+1):
        ls.append(list())
    ls[0].append(dickson(n,begin))
    cur = begin
    for i in range(1,n+1):
        cur+=h
        ls[i].append(dickson(n,cur))
    for j in range(1,n+1):
        for k in range(n-j,-1,-1):
            ls[k].append(ls[k+1][j-1]-ls[k][j-1])
    return ls[0]


'''
    Schonhage-Strassen

From now on we define functions to implement the
schonhage-strassen polynomial multiplication algorithm in the ring R = Z/nZ
'''

'''
Multiply a polynomial P by the power of a principal root of unity x^k
Almost like a cyclic shift
'''
def mulX(P,k,M):
    if P == []: return [0]*M
    if len(P)<M:
        P+=[0]*(M-len(P))
    t,kb = k//M, k%M
    if t%2==0:
        return [-x for x in P[M-kb:M]]+P[:M-kb]
    else:
        return P[M-kb:M]+[-x for x in P[:M-kb]]

'''
Returns P mod X^m+1
'''
def remn1(P,m,n):
    res = P[:m]
    nex = P[m:2*m]
    ind = m
    minu = True
    while nex!=[]:
        if minu:
            res = minus(res,nex,n)
        else:
            res = plus(res,nex,n)
        minu = not minu
        ind+=m
        nex = P[ind:ind+m]
    return res

'''
add P and Q
'''
def plus(P,Q,n):
    if len(P)>len(Q):
        l = len(P)
        res = [0]*l
        for i in range(len(Q)):
            res[i] = (P[i]+Q[i])%n
        for i in range(len(Q),l):
            res[i] = P[i]
    else:
        l = len(Q)
        res = [0]*l
        for i in range(len(P)):
            res[i] = (P[i]+Q[i])%n
        for i in range(len(P),l):
            res[i] = Q[i]
    return res

'''
Discrete Fourier transform in the ring
D[y] where D = R[x]/(x^M+1)
'''
def DFT(P,k,M,u,n):
    if u==1:
        return [P[0]]
    t = u>>1
    r0 = [[]]*t
    for j in range(t):
        r0[j] = plus(P[j],P[j+t],n)
    r1 = [[]]*t
    i = 0
    for j in range(t):
        r1[j] = mulX(minus(P[j],P[j+t],n),i,M)
        i += k
    res0 = DFT(r0,k<<1,M,t,n)
    res1 = DFT(r1,k<<1,M,t,n)
    res = [0]*u
    res[::2] =  res0
    res[1::2] = res1
    return res

def fast_convolution(P,Q,k,w,M,m,n):
    a = DFT(P,w,M,k,n)
    b = DFT(Q,w,M,k,n)
    res = [[]]*k
    for i in range(k):
        res[i] = schon_strass_recur(a[i], b[i],m,n)
    inv2 = (n+1)>>1
    invk = powmod(inv2,int(log(k,2)),n)
    return [[(invk*c)%n for c in p] for p in DFT(res,2*M-w,M,k,n)]

def mul(P,Q,n):
    resLen = len(P)+len(Q)-1
    recurLen = (resLen-1).bit_length()
    return schon_strass_recur(P, Q, recurLen,n)[:resLen]


def schon_strass_recur(P,Q,s,n):
    if s<10:
        return remn1(mul_kronecker(P,Q,n),1<<s,n)
    if s%2==0:
        m = t = s//2
        eta = 2
    else:
        m = s//2
        t = m+1
        eta = 1
    
    s1,m1,t1 = 1<<s, 1<<m, 1<<t
    P_ = []
    Q_ = []
    j = 0
    for i in range(0,s1,m1):
        P_.append(mulX(P[i:i+m1],j,2*m1))
        Q_.append(mulX(Q[i:i+m1],j,2*m1))
        j+=eta
    res = fast_convolution(P_, Q_, t1, 2*eta, 2*m1,m+1,n)
    j = eta
    for i in range(1,len(res)):
        res[i] = mulX(res[i],4*m1-j,2*m1)
        j+=eta

    res2 = res[0][:m1]
    z = m1
    for i in range(len(res)-1):
        res2[z+1:z+m1+1] = plus(res[i][m1:], res[i+1][:m1], n)
        z += m1
    res2[z:z+m1] = res[-1][m1:]
    return remn1(res2,s1,n)

class CachedFFT():
    def __init__(self,P,s,n):
        self.s = s
        self.s1 = 1<<s
        if s<10:
            self.P = P
            return
        if s%2==0:
            self.m = self.t = s//2
            self.eta = 2
        else:
            self.m = s//2
            self.t = self.m+1
            self.eta = 1
        self.m1,self.t1 = 1<<self.m, 1<<self.t
        P_ = []
        j = 0
        for i in range(0,self.s1,self.m1):
            P_.append(mulX(P[i:i+self.m1],j,2*self.m1))
            j+=self.eta
        self.transforms = []
        for poly in DFT(P_, 2*self.eta, 2*self.m1, self.t1, n):
            self.transforms.append(CachedFFT(poly,self.m+1,n))


    