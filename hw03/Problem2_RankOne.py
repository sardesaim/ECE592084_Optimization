# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 10:59:51 2019

@author: sarde
"""
#Problem 1 
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d

def Bracketing(f, alphInit, epsi):
    x0 = alphInit
    x1 = x0+epsi
    x2 = x0+ 2*epsi
    n = r = 1
    if(f(x0)>f(x1) and f(x1)<f(x2)):
        a0 = x0
        b0 = x2
    else:
        while not (f(x0)>f(x1) and f(x1)<f(x2)):
            if((f(x0)>f(x1) and f(x1)>f(x2))):
                x0 = x1
                x1 = x2
                x2 = x2 +(2**(n+1)) * epsi
                n+=1
            elif((f(x0)<f(x1) and f(x1)<f(x2))):
                x2 = x1
                x1 = x0
                x0 = x0 - (2**r) * epsi
                r+=1
        a0 = x0 
        b0 = x2
        return a0,b0

def goldSec(a0,b0,fcn):
#    epsi = 0.075
    N = np.ceil(np.log(0.001/np.linalg.norm(b0-a0))/np.log(.6180))
    rho = .382
    a = a0
    b = b0
    s = a+rho*(b-a)
    t = a + (1-rho) * (b-a)
    f1 = fcn(s)
    f2 = fcn(t)
    for i in range(1,int(N+1)):
        if(fcn(s)<fcn(t)):
            b = t
            t = s
            s = a+rho*(b-a)
            f2 = f1
            f1 = fcn(s)
        elif(fcn(s)>fcn(t)):
            a = s
            s = t
            t = a+(1-rho)*(b-a)
            f1 = f2
            f2 = fcn(t)
        else:
            break
    return s,t

def RankOne(X_init, f):
    normx = 1
    H = []
    H.append(np.array([[10,2],[2,10]]))
    if grad(X_init[0][0], X_init[0][1]) == (0,0):
        return X_init
    else:        
        Hk = H[-1]
        gradk = grad(X_init[0][0],X_init[0][1])
        d = -Hk.dot(gradk)
        fAlph = lambda alpha: (f(X_init[-1][0] + alpha*d[0], X_init[-1][0] + alpha*d[1]))
        al1, al2 = Bracketing(fAlph, 0.2, 0.05)
        alpa, alpb = goldSec(al1, al2, fAlph)
        alpha = (alpa+alpb)/2
        X_init.append(X_init[-1]+np.multiply(alpha,d))
        iters = 1
        while normx > epsilon and iters < max_iterations:
            gradk = grad(X_init[-1][0],X_init[-1][1])
            dk = -Hk.dot(gradk)
            fAlph = lambda alpha: (f(X_init[-1][0] + alpha*dk[0], X_init[-1][1] + alpha*dk[1]))
            al1, al2 = Bracketing(fAlph, 0.2, 0.05)
            alpa, alpb = goldSec(al1, al2, fAlph)
            alphak = (alpa+alpb)/2
            X_init.append(X_init[iters]+np.multiply(alphak,dk))
            deltaXk = alphak*dk
            deltaGk = np.subtract(grad(X_init[-1][0],X_init[-1][1]),gradk)
            H.append(np.array(Hk+(((deltaXk-Hk.dot(deltaGk)).dot(np.transpose((deltaXk-Hk.dot(deltaGk)))) )\
                         /(np.transpose(deltaGk).dot(deltaXk-Hk.dot(deltaGk))))))
            normx = np.linalg.norm(X_init[-1]-X_init[-2])
            iters+=1
        return X_init
    
def DFP(X_init, f):
    normx = 1
    H = []
    H.append(np.array([[10,2],[2,10]]))
    if grad(X_init[0][0], X_init[0][1]) == (0,0):
        return X_init
    else:        
        Hk = H[-1]
        gradk = grad(X_init[0][0],X_init[0][1])
        d = -Hk.dot(gradk)
        fAlph = lambda alpha: (f(X_init[-1][0] + alpha*d[0], X_init[-1][0] + alpha*d[1]))
        al1, al2 = Bracketing(fAlph, 0.2, 0.05)
        alpa, alpb = goldSec(al1, al2, fAlph)
        alpha = (alpa+alpb)/2
        X_init.append(X_init[-1]+np.multiply(alpha,d))
        iters = 1
        while normx > epsilon and iters < max_iterations:
            gradk = grad(X_init[-1][0],X_init[-1][1])
            dk = -Hk.dot(gradk)
            fAlph = lambda alpha: (f(X_init[-1][0] + alpha*dk[0], X_init[-1][1] + alpha*dk[1]))
            al1, al2 = Bracketing(fAlph, 0.2, 0.05)
            alpa, alpb = goldSec(al1, al2, fAlph)
            alphak = (alpa+alpb)/2
            X_init.append(X_init[iters]+np.multiply(alphak,dk))
            deltaXk = alphak*dk
            deltaGk = np.subtract(grad(X_init[-1][0],X_init[-1][1]),gradk)
            H.append(np.array(Hk+(deltaXk.dot(np.transpose(deltaXk))/(np.transpose(deltaXk).dot(deltaGk))) -\
                              (np.transpose((Hk.dot(deltaGk))).dot(Hk.dot(deltaGk)))/ \
                              ((np.transpose(deltaGk).dot(Hk)).dot(deltaGk))))
            normx = np.linalg.norm(X_init[-1]-X_init[-2])
            iters+=1
        return X_init

def BFGS(X_init, f):
    normx = 1
    H = []
    H.append(np.array([[10,2],[2,10]]))
    if grad(X_init[0][0], X_init[0][1]) == (0,0):
        return X_init
    else:        
        Hk = H[-1]
        gradk = grad(X_init[0][0],X_init[0][1])
        d = -Hk.dot(gradk)
        fAlph = lambda alpha: (f(X_init[-1][0] + alpha*d[0], X_init[-1][0] + alpha*d[1]))
        al1, al2 = Bracketing(fAlph, 0.2, 0.05)
        alpa, alpb = goldSec(al1, al2, fAlph)
        alpha = (alpa+alpb)/2
        X_init.append(X_init[-1]+np.multiply(alpha,d))
        iters = 1
        while normx > epsilon and iters < max_iterations:
            gradk = grad(X_init[-1][0],X_init[-1][1])
            dk = -Hk.dot(gradk)
            fAlph = lambda alpha: (f(X_init[-1][0] + alpha*dk[0], X_init[-1][1] + alpha*dk[1]))
            al1, al2 = Bracketing(fAlph, 0.2, 0.05)
            alpa, alpb = goldSec(al1, al2, fAlph)
            alphak = (alpa+alpb)/2
            X_init.append(X_init[iters]+np.multiply(alphak,dk))
            deltaXk = alphak*dk
            deltaGk = np.subtract(grad(X_init[-1][0],X_init[-1][1]),gradk)
            H.append(np.array(Hk+(((deltaXk-Hk.dot(deltaGk)).dot(np.transpose((deltaXk-Hk.dot(deltaGk)))) )\
                         /(np.transpose(deltaGk).dot(deltaXk-Hk.dot(deltaGk))))))
            normx = np.linalg.norm(X_init[-1]-X_init[-2])
            iters+=1
        return X_init

def plotSeq(X_init):
    yplt = xplt = np.arange(-1, 1, 0.025)
    Xpt, Ypt = np.meshgrid(xplt, yplt)
    Z = f(Xpt,Ypt)
    fig, ax = plt.subplots()
    CS = ax.contour(Xpt, Ypt, Z)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title('Sequence of points')
    plx = []
    ply = []
    for i in range(len(X_init)):
        plx.append(X_init[i][0])
        ply.append(X_init[i][1])
    plt.plot(plx, ply , 'b-^')

if __name__ == "__main__":
    #define function in terms of x1 x2 
    f = lambda x1, x2: ((x2-x1)**4 + 12*x1*x2 - x1 + x2 - 3)
    #supply gradients 
    grad = lambda x1, x2: ([(12*x2 + 4*(x1-x2)**3 - 1),(12*x1 - 4*(x1-x2)**3 + 1)])
    #dy = lambda x1, x2: 
    alpha = 0.05
    epsilon = 0.001
    max_iterations = 1000
    X_init = []
    X_init.append(np.array([0.55,0.7]))
    X_init = RankOne(X_init,f)
    plotSeq(X_init)
    
    X_init1 = []
    X_init1.append(np.array([-0.9,-0.5]))
    X_init1 = RankOne(X_init1,f)
    plotSeq(X_init1)
    
    X_init2 = []
    X_init2.append(np.array([-0.9,-0.5]))
    X_init2 = DFP(X_init2,f)
    plotSeq(X_init2)