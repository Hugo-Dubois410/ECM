# -*-coding:UTF-8 -*

'''
Created on 5 févr. 2017
@author: garzol
'''

from factorization import factor, sieve
import sys, time
from multiprocessing import freeze_support
freeze_support()

def get_arg(n):
    try: 
        n = int(n)
    except: 
        try:
            n = eval(n)
            if(type(n).__name__ != "int"):
                raise ValueError
        except:
            raise ValueError
    return n

def main():
    print("\nECM Factoring method using python3 and gmpy2")
    PRIMES = sieve(120000000)
    nProc = 1
    B1 = 100000
    B2 = 450000000
    print("B1 =",B1)
    print("B2 =",B2)
    print("Number of parallel curves:",nProc)
    print("\nEnter an integer or an integer valued expression to try and factor it.")
    print("You can also change the parameters by typing param=value\nWhere param is one of B1, B2 or curves")
    while True:
        sys.stdout.write("\n>> ")
        n = input()
        if n=="quit":
            sys.exit(0)
        elif n.startswith("B1="):
            try:
                B1 = get_arg(n[3:])
            except ValueError: continue
            print("\nB1 =",B1)
            print("B2 =",B2)
            print("Number of curves:",nProc)
        elif n.startswith("B2="):
            try:
                B2 = get_arg(n[3:])
            except ValueError: continue
            print("\nB1 =",B1)
            print("B2 =",B2)
            print("Number of curves:",nProc)
        elif n.startswith("curves="):
            try:
                nProc = get_arg(n[7:])
            except ValueError: continue
            print("\nB1 =",B1)
            print("B2 =",B2)
            print("Number of curves:",nProc)
        else:
            try: 
                n = get_arg(n)
            except ValueError:
                print(n,"is not an integer")
                continue
            startt = time.time()
            print(factor(n,PRIMES,nProc,B1,B2),"are the (probable) prime factors of",n)
            print("elapsed time",time.time()-startt)

if __name__=='__main__':
    try : main()
    except (EOFError, KeyboardInterrupt):
        sys.exit(0)
    

