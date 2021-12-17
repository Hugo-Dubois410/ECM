# -*-coding:UTF-8 -*

from gmpy2 import gcdext

def invert(xs,n):
    y = xs[0]
    cs = [y]
    for x in xs[1:]:
        y=(y*x)%n
        cs.append(y)
    g, inv, _ = gcdext(y,n)
    if g!=1 : 
        return g
    cs[-1] = (cs[-2]*inv)%n
    for i in range(len(cs)-2,0,-1):
        inv = (inv*xs[i+1])%n
        cs[i] = (cs[i-1]*inv)%n
    cs[0] = (inv*xs[1])%n
    return cs
    


class EllCurveEdwards(object):

    def __init__(self,d,n):
        self.d, self.d2, self.n = d , d<<1, n 
                 
    def add_unified(self,P,Q):
        [X1,Y1,Z1,T1] = P
        [X2,Y2,Z2,T2] = Q
        A = ((Y1-X1)*(Y2-X2))%self.n
        B = ((Y1+X1)*(Y2+X2))%self.n
        C = (((T1*self.d2)%self.n)*T2)%self.n
        D = (Z1*2*Z2)%self.n
        E = (B-A)%self.n
        F = (D-C)%self.n
        G = (D+C)%self.n
        H = (B+A)%self.n
        X3 = (E*F)%self.n
        Y3 = (G*H)%self.n
        T3 = (E*H)%self.n
        Z3 = (F*G)%self.n
        return [X3,Y3,Z3,T3]
    
    def add(self,P,Q):
        [X1,Y1,Z1,T1] = P
        [X2,Y2,Z2,T2] = Q
        A = ((Y1-X1)*(Y2+X2))%self.n
        B = ((Y1+X1)*(Y2-X2))%self.n
        C = (Z1*2*T2)%self.n
        D = (T1*2*Z2)%self.n
        E = (D+C)%self.n
        F = (B-A)%self.n
        G = (B+A)%self.n
        H = (D-C)%self.n
        X3 = (E*F)%self.n
        Y3 = (G*H)%self.n
        T3 = (E*H)%self.n
        Z3 = (F*G)%self.n
        return [X3,Y3,Z3,T3]

    def double(self,P):
        [X1,Y1,Z1,_] = P
        A = (X1**2)%self.n
        B = (Y1**2)%self.n
        C = (2*(Z1**2))%self.n
        D = -A
        E = (((X1+Y1)**2)-A-B)%self.n
        G = (D+B)%self.n
        F = (G-C)%self.n
        H = (D-B)%self.n
        X3 = (E*F)%self.n
        Y3 = (G*H)%self.n
        T3 = (E*H)%self.n
        Z3 = (F*G)%self.n
        return [X3,Y3,Z3,T3]
    
    def mul(self,k,P):
    #binary multiplication 
        result, power_ = [0,1,1,0], P
        while k>0:
            if k & 1:
                result = self.add(result,power_)
            k>>=1
            power_ = self.double(power_)
        return result
    
    def toWeierstrass(self,P):
        g, z ,_ = gcdext(P[2],self.n)
        if g!=1:
            return g
        x = (P[0]*z)%self.n
        y = (P[1]*z)%self.n
        g2, t ,_ = gcdext(self.d+1,self.n)
        if g2!=1:
            return g2
        A = (2*(-self.d+1)*t)%self.n
        B = (-4*t)%self.n
        g3, y1 ,_ = gcdext(1-y,self.n)
        if g3!=1:
            return g3
        g4, x1 ,_ = gcdext(x,self.n)
        if g4!=1:
            return g4
        u = ((1+y)*y1)%self.n
        v = (u*x1)%self.n
        g5, invB ,_ = gcdext(B,self.n)
        if g5!=1:
            return g5
        _, inv3 ,_ = gcdext(3,self.n)
        _,  inv27 ,_ = gcdext(27,self.n)
        i, j = ((u*invB)+A*inv3*invB)%self.n, (v*invB)%self.n
        a, b = ((3-A*A)*inv3*invB*invB)%self.n, ((2*A*A*A-9*A)*invB*invB*invB*inv27)%self.n
        return [EllCurve(a,b,self.n), (i,j,1)]
    
class EllCurve(object):

    def __init__(self,a, b ,n):
        self.a , self.b, self.n =a, b , n 
                 
    def add(self,P,Q):
    
        if P[2]==0:
            return Q
        if Q[2]==0:
            return P
        if P[0]==Q[0]:
            if P[1] != Q[1]:
                return [0,0,0]
            d = P[1]<<1
            pgcd, inv,_ = gcdext(d,self.n)
            if pgcd !=1 : return [0,0,pgcd]
            x = (3*((P[0]*P[0])%self.n)+self.a)%self.n
            slope = (x*inv)%self.n
            x = (slope*slope-(P[0]<<1))%self.n
            return [x,(slope*(P[0]-x)-P[1])%self.n, 1 ]
        d = P[0]-Q[0]
        pgcd, inv, _  = gcdext(d,self.n)
        if pgcd !=1 : return [0,0,pgcd]
        slope = ((P[1]-Q[1])*inv)%self.n
        x = (slope*slope-P[0]-Q[0])%self.n
        return [ x, (slope*(x-P[0])+P[1])%self.n,1]
    

    def add_parallel(self,P):
        ds = list()
        for Q,R in zip(P,P[1:]):
            ds.append(Q[0]-R[0])
        inverted = invert(ds,self.n)
        if not isinstance(inverted, list):
            return inverted
        for i in range(len(P)-1):
            slope = ((P[i][1]-P[i+1][1])*inverted[i])%self.n
            P[i][0] = (slope * slope - P[i][0] - P[i+1][0]) % self.n
            P[i][1] = (slope * (P[i][0] - P[i+1][0]) + P[i+1][1]) % self.n
    
    def mul(self,k,P):
    #binary multiplication 
        result, power_ = [0,0,0], P
        while k>0:
            if k & 1:
                result = self.add(result,power_)
                if result[2]>1:
                    return result
            k>>=1
            power_ = self.add(power_, power_)
        return result
    
