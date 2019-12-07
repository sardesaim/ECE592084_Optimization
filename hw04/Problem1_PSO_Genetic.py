# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 10:59:51 2019

@author: sarde
"""
#Problem 1 
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d

class GlobalSearchAlgos():
    #define function in terms of x1 x2 
    epsilon = 0.02
    max_iterations = 1000
    numPositions = 0
    epochs = 10
    gbvals=[]
    wrsVals = []
    meanVals = []
    def __init__(self,d):
        self.numPositions = d
    def f(self,x1,x2):
        f = ((x2-x1)**4 + 12*x1*x2 - x1 + x2 - 3)
        return f
    def PSO(self):
        """
        Position Swarm Optimization
        inputs - 
        outputs - gbest
        """
        w = 0.9
        c1 = 2
        c2 = 2
        X_pos = []
        p_pos = []
        vel = []
        X_pos.append(np.random.rand(self.numPositions,2)*2.0-1)
        vel.append(np.random.rand(self.numPositions,2))
        k = 0
        fval = []
        for i,val in enumerate(X_pos[-1]):
            fval.append(self.f(X_pos[-1][i][0],X_pos[-1][i][1]))
        gBest = []
        p_pos = X_pos
        gBest.append(X_pos[-1][fval.index(min(fval))])
        while k<self.epochs:
            r=np.random.rand(self.numPositions,2)
            s=np.random.rand(self.numPositions,2)    
            vel.append(np.multiply(w,vel[-1])+np.multiply(np.multiply(c1,r), (p_pos[-1]-X_pos[-1])) + \
                       np.multiply(np.multiply(c2,s), (gBest[-1]-X_pos[-1])))
            p_pos.append(p_pos[-1])
            X_pos.append(X_pos[-1] + vel[-1])
            fval=[]
            gBestTemp = gBest[-1]
            for i,val in enumerate(X_pos[-1]):
                fval.append(self.f(X_pos[-1][i][0],X_pos[-1][i][1]))
                if(fval[-1]<self.f(p_pos[-1][i][0],p_pos[-1][i][1])):
                    p_pos[-1][i] = X_pos[-1][i]
                else:
                    pass
                if(fval[-1]<self.f(gBest[-1][0],gBest[-1][1])):
                    gBestTemp = X_pos[-1][i]
                else:
                    pass
            gBest.append(gBestTemp)
            self.gbvals.append(self.f(gBest[-1][0], gBest[-1][1]))
            self.meanVals.append(np.mean(fval))
            self.wrsVals.append(self.f(X_pos[-1][fval.index(max(fval))][0],X_pos[-1][fval.index(max(fval))][1]))
#            print(fval,end='\n')
            k+=1
        return gBest[-1]
    
    def bit2num(self,bit, rng):
        integ = np.polyval(bit,2)
        num = rng[0]+integ*(rng[1]-rng[0])/(2^(len(bit)-1))
        return num
    
    def evalPopu(self, popu, bit_n, rng, obj_func):
        popu_s = np.shape(popu)[0]
        str_len = np.shape(popu)[1]
        fitness = np.zeros(popu_s)
        for i in range(popu_s):
            num1 = self.bit2num(popu[i][0:bit_n], rng)
            num2 = self.bit2num(popu[i][bit_n+1:str_len], rng)
            fitness[i] = obj_func(num1, num2)
        return fitness
    
    def nextPopu(self,popu, fitness, x_over, mut_rate):
        new_popu = popu 
        popu_s = np.shape(popu)[0]
        str_len = np.shape(popu)[1]
        #rescaling fitness
        fitness = fitness - min(fitness)
        tmp_fitness = fitness
        index1 = np.argmax(tmp_fitness)
        tmp_fitness[index1]=tmp_fitness[np.argmin(tmp_fitness)]
        index2 = np.argmax(tmp_fitness)
        new_popu[0:1][:] = popu[index1:index2][:]
        #for roulette wheel
        total = sum(fitness)
        if total == 0:
            fitness = np.ones(popu_s,1)/popu_s
        cumprob = np.cumsum(fitness)/total
        #selection and crossover
        rnd = np.random.rand()
        for i in range(2, popu_s/2):
            tmp = np.argwhere(cumprob-rnd>0)
            parent1 = popu[tmp[1]][:]
            tmp = np.argwhere(cumprob-rnd>0)
            parent2 = popu[tmp[1]][:]
            #do crossover
            if rnd>x_over:
                x_over_pt = np.ceil(rnd*(str_len-1))
                new_popu[i*2-1][:] = np.hstack(parent1[1:x_over_pt],parent2[x_over_pt+1:str_len-1])
                new_popu[i*2][:] = np.hstack(parent2[1:x_over_pt],parent1[x_over_pt+1:str_len-1])
        #mutation elites are not subject to this 
        mask = np.random.rand(popu_s, str_len)<mut_rate
        new_popu = new_popu ^ mask
        return new_popu
    
    def CanonicalGA(self):
        """
        Position Swarm Optimization
        inputs - 
        outputs - gbest
        """
        gen_n = 50; 
        popu_size = 20; 
        xover_rate = 0.85
        mutate_rate = 0.01
        bit_n = 16 
        obj_func = lambda x1, x2: self.f(x1, x2)
        var_n = 2
        rng = np.array([-1, 1])
        popu = np.random.rand(popu_size, bit_n*var_n)>0.5;
        upper = np.zeros((gen_n, 1))
        average = np.zeros((gen_n, 1))
        lower  = np.zeros((gen_n, 1))
        for i in range(gen_n):
            func_val = self.evalPopu(popu, bit_n, rng, obj_func)
            upper[i] = np.max(func_val)
            average[i] = np.mean(func_val)
            lower[i] = np.min(func_val)
            print(lower[i],end='\n')
            popu = self.nextPopu(popu, func_val, xover_rate, mutate_rate)
            
    def plotCont(self):
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
        ax.set_title('Global Minima on Contour')
    
    def plotObjVal(self):
        """
        Plot the sequence of Xk on level sets of f 
        and save the level sets and points in a .csv file
        """
        plt.figure()
        x = np.linspace(1,self.epochs, self.epochs)
        plt.plot(x,self.gbvals, 'r', x, self.wrsVals, 'b', x, self.meanVals, 'g')
#        plt.plot(x,self.wrsVals)
#        plx = []
#        ply = []
#        fval = []
#        for i in range(len(X_init)):
#            plx.append(X_init[i][0])
#            ply.append(X_init[i][1])
#            fval.append(self.f(plx[-1], ply[-1]))
#        plt.plot(plx, ply , 'b-^')
        
if __name__ == "__main__":  
    ps = GlobalSearchAlgos(20)
    gB = ps.PSO()
    plt.figure()
    ps.plotCont()
    plt.plot(gB[0], gB[1] , 'bo-')
    ps.plotObjVal()
#    ps.CanonicalGA()