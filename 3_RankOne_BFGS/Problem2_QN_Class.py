# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 10:59:51 2019

@author: sarde
"""
#Problem 1 
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d

class QuasiNewton():
    #define function in terms of x1 x2 
    epsilon = 0.02
    max_iterations = 1000
    X_init = []
    def __init__(self, X_init):
        self.X_init=X_init
    def f(self,x1,x2):
        f = ((x2-x1)**4 + 12*x1*x2 - x1 + x2 - 3)
        return f
    def grad(self,x1,x2):
        grad = ([(12*x2 + 4*(x1-x2)**3 - 1),(12*x1 - 4*(x1-x2)**3 + 1)])
        return grad
    def Bracketing(self, f, alphInit, epsi):
        """
        Does bracketing to set up initial interval for line search
        inputs - function, initial parameter, epsilon
        outputs - initial interval 
        """
        x0 = alphInit
        x1 = x0+epsi
        x2 = x0+ 2*epsi
        n = r = 1
        if(f(x0)>f(x1) and f(x1)<f(x2)):
            a0 = x0
            b0 = x2
            return a0,b0
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
    
    def goldSec(self, a0,b0,fcn):
        """
        Reduces uncertainty range of the interval found out by bracketing 
        alpha can be picked from the reduced uncertainty interval
        inputs - bracketed interval, function 
        outputs - reduced uncertainty interval
        """
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
        return a,b
    
    def RankOne(self, X_init, f):
        """
        Basic rank one update for estimating the Inverse hessian 
        Hk 
        Initial Hk chosen to be 10,2; 2, 10 - pd and symmetric 
        inputs - Initial array of X, function f
        outputs - updated array of X 
        """
        normx = 1
        H = []
        H.append(np.array([[1,0],[0,1]]))
        #if gradient of initial point is zero STOP
        if self.grad(X_init[-1][0], X_init[-1][1]) == (0,0):
            return X_init
        else:        
            Hk = H[-1]  #set Hk as last element in H[]
            gradk = self.grad(X_init[0][0],X_init[0][1]) #find kth gradient
            d = -Hk.dot(gradk)  #set direction according to kth gradient and Hk
            fAlph = lambda alpha: (f(X_init[-1][0] + alpha*d[0], X_init[-1][0] \
                                     + alpha*d[1])) #find a fAlph 
            al1, al2 = self.Bracketing(fAlph, 0, 1) #find optimum alpha interval
            alpa, alpb = self.goldSec(al1, al2, fAlph) 
            alphak = (alpa+alpb)/2  #find optimum alpha 
            X_init.append(X_init[-1]+np.multiply(alphak,d)) #find second point
            iters = 1
            #iterate till consecutive points are within a small value epsilon
            while normx > self.epsilon and self.grad(X_init[-1][0],X_init[-1][1]) and iters < self.max_iterations:   
                #update Hk as the last element in the H estimates list
                Hk = H[-1]
                gradk = self.grad(X_init[-1][0],X_init[-1][1]) #find gradient of latest pt
                dk = -Hk.dot(gradk) #find the direction from Hk and gk
                #line search over x_alpha*d
                fAlph = lambda alphak: (f(X_init[-1][0] + alphak*dk[0], \
                                          X_init[-1][1] + alphak*dk[1]))    
                #bracketing and golden section to find optimum alpha
                al1, al2 = self.Bracketing(fAlph, 0, 1)
                alpa, alpb = self.goldSec(al1, al2, fAlph)
                alphak = (alpa+alpb)/2
                #append updated point to the list of points
                X_init.append(X_init[iters]+np.multiply(alphak,dk))
                deltaXk = alphak*dk 
                deltaGk = np.subtract(self.grad(X_init[-1][0],X_init[-1][1]),gradk)
                #Rank one update of Hk
                H.append(np.array(Hk+(((deltaXk-Hk.dot(deltaGk)).dot\
                                       (np.transpose((deltaXk-Hk.dot(deltaGk)))) )\
                             /(np.transpose(deltaGk).dot(deltaXk-Hk.dot(deltaGk))))))
                normx = np.linalg.norm(X_init[-1]-X_init[-2])
                #if np.all(np.linalg.eigvals(H[-1]) > 0):
                #    print('The estimated Hk+1 is not pd, descent not guaranteed')
                #    break
                #else:
                #    continue
                iters+=1
            return X_init
        
    def plotSeq(self, X_init, meth_eg):
        """
        Plot the sequence of Xk on level sets of f 
        and save the level sets and points in a .csv file
        """
        yplt = xplt = np.arange(-1, 1, 0.025)
        Xpt, Ypt = np.meshgrid(xplt, yplt)
        Z = self.f(Xpt,Ypt)
        fig, ax = plt.subplots(num = None, figsize = (8,6), dpi = 90, facecolor = 'w', edgecolor = 'k')
    #    plt.figure()
        CS = ax.contour(Xpt, Ypt, Z)
        ax.clabel(CS, inline=1, fontsize=10)
        ax.set_title('Sequence of points')
        plx = []
        ply = []
        fval = []
        for i in range(len(X_init)):
            plx.append(X_init[i][0])
            ply.append(X_init[i][1])
            fval.append(self.f(plx[-1], ply[-1]))
        plt.plot(plx, ply , 'b-^')
        if meth_eg == 1:
            np.savetxt('RankOneSeqPtsEg1.csv', np.column_stack((plx,ply,fval)), \
                       delimiter = ",", fmt = '%s')
        elif meth_eg == 2:
            np.savetxt('RankOneSeqPtsEg2.csv', np.column_stack((plx,ply,fval)), \
                       delimiter = ",", fmt = '%s')
        elif meth_eg == 3:
            np.savetxt('DFPSeqPtsEg1.csv', np.column_stack((plx,ply,fval)), \
                       delimiter = ",", fmt = '%s')
        elif meth_eg == 4:
            np.savetxt('DFPSeqPtsEg2.csv', np.column_stack((plx,ply,fval)), \
                       delimiter = ",", fmt = '%s')
        elif meth_eg == 5:
            np.savetxt('BFGSSeqPtsEg1.csv', np.column_stack((plx,ply,fval)), \
                       delimiter = ",", fmt = '%s')
        elif meth_eg == 6:
            np.savetxt('BFGSSeqPtsEg2.csv', np.column_stack((plx,ply,fval)), \
                       delimiter = ",", fmt = '%s')
    #    plt.quiver(plx[:-1], ply[:-1], scale_units='xy', angles='xy', scale=1)

if __name__ == "__main__":  
    X = []
    X.append(np.array([0.55,0.7]))
    qn = QuasiNewton(X)
    qn.X_init = qn.RankOne(qn.X_init,qn.f)
    qn.plotSeq(qn.X_init,1)
    