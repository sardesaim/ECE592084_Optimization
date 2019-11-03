# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 00:31:54 2019

@author: sarde
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#define function in terms of x1 x2 
f = lambda x1, x2: (x1**2 + x1*x2 + x2**2)
X = []
X.append(np.array([0.8,-0.25]))
epsi = 0.075
k = 1

def Bracketing(f, X, epsi):
    x0 = X[0]
    x1 = x0+epsi
    x2 = x0+ 2*epsi
    n = r = 1
    if(f(x0[0], x0[1])>f(x1[0], x1[1]) and f(x1[0], x1[1])<f(x2[0], x2[1])):
        a0 = x0
        b0 = x2
    else:
        while not (f(x0[0], x0[1])>f(x1[0], x1[1]) and f(x1[0], x1[1])<f(x2[0], x2[1])):
            if((f(x0[0], x0[1])>f(x1[0], x1[1]) and f(x1[0], x1[1])>f(x2[0], x2[1]))):
                x0 = x1
                x1 = x2
                x2 = x2 +(2**(n+1)) * epsi
                n+=1
            elif((f(x0[0], x0[1])<f(x1[0], x1[1]) and f(x1[0], x1[1])<f(x2[0], x2[1]))):
                x2 = x1
                x1 = x0
                x0 = x0 - (2**r) * epsi
                r+=1
        a0 = x0 
        b0 = x2
        return a0,b0

def goldSec(a0,b0,fcn):
#    epsi = 0.075
    N = np.ceil(np.log(0.01/np.linalg.norm(b0-a0))/np.log(.6180))
    rho = .382
    a = a0
    b = b0
    s = a+rho*(b-a)
    t = a + (1-rho) * (b-a)
    f1 = fcn(s[0], s[1])
    f2 = fcn(t[0], t[1])
    for i in range(1,int(N+1)):
        if(fcn(s[0], s[1])<fcn(t[0], t[1])):
            b = t
            t = s
            s = a+rho*(b-a)
            f2 = f1
            f1 = fcn(s[0], s[1])
        elif(fcn(s[0], s[1])>fcn(t[0], t[1])):
            a = s
            s = t
            t = a+(1-rho)*(b-a)
            f1 = f2
            f2 = fcn(t[0], t[1])
        else:
            break
    return s,t
    
a0,b0 = Bracketing(f, X, epsi)
(s,t) = goldSec(a0,b0,f)