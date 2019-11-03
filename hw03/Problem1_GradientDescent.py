# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 10:59:51 2019

@author: sarde
"""
#Problem 1 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#define function in terms of x1 x2 
f = lambda x1, x2: ((x2-x1)**4 + 12*x1*x2 - x1 + x2 - 3)
#supply gradients 
grad = lambda x1, x2: ((12*x2 + 4*(x1-x2)**3 - 1),(12*x1 - 4*(x1-x2)**3 + 1))
#dy = lambda x1, x2: 

alpha = 0.02
epsilon = 0.001
max_iterations = 1000
normx = 1
X_init = []
X_init.append(np.array([0.55,0.7]))
X_init.append(X_init[-1]-np.multiply(alpha,grad(X_init[0][0], X_init[0][1])))
iters = 1
while normx > epsilon and iters < max_iterations:
    X_init.append(X_init[iters]-np.multiply(alpha,grad(X_init[iters][0], X_init[iters][1])))
    normx = np.linalg.norm(X_init[-1]-X_init[-2])
    iters+=1
    
#yplt = xplt = np.linspace(-1, 1, 100)
#np.meshgrid(xplt, yplt)
#fig = plt.figure()
#ax = plt.axes(projection = '3d')
#ax.plot3D(xplt, yplt, f(xplt, yplt))