#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 00:30:20 2022

@author: nrobinson15
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt



G= 6.67408*(10**-11)
M= 1.989*(10**30)

x1 = 2.52*1.496*10**11
y1 = 0
v_x1 = 0
v_y1 = sqrt(G*M/x1)

x2 = 5.24*1.496*10**11
y2 = 0
v_x2 = 0
v_y2 = sqrt(G*M/x2)


x3 = 0
y3 = 0
v_x3 = 0
v_y3 = 0

m1 = (10**-3)*M
m2 = (4*(10**-2))*M
m3 = M

k=[] 

# Add extra terms on original functions, add 4 more functions for the specific planet, 
# change the intial conditions with a1 = 2.52 AU and a2 = 5.24 AU, How to find initial velocities of planets? Look it up



def first_order(k, t):
    Fx12 = G*m1*m2*(k[4]-k[0]) / (((k[4]-k[0])**2 + (k[5]-k[1])**2)**1.5)
    Fx21 = -Fx12
    
    Fx31 = G*m1*m3*(k[0]-k[8]) / (((k[0]-k[8])**2 + (k[1]-k[9])**2)**1.5)
    Fx13 = -Fx31
    
    Fx32 = G*m3*m2*(k[4]-k[8]) / (((k[4]-k[8])**2 + (k[5]-k[9])**2)**1.5)
    Fx23 = -Fx32
    
    Fy12 = G*m1*m2*(k[5]-k[1]) / (((k[4]-k[0])**2 + (k[5]-k[1])**2)**1.5)
    Fy21 = -Fy12
    
    Fy31 = G*m1*m3*(k[1]-k[9]) / (((k[0]-k[8])**2 + (k[1]-k[9])**2)**1.5)
    Fy13 = -Fy31
    
    Fy32 = G*m3*m2*(k[5]-k[9]) / (((k[4]-k[8])**2 + (k[5]-k[9])**2)**1.5)
    Fy23 = -Fy32
    
    
    
    return np.array([k[2],
                     k[3],
                     (Fx12+Fx13)/m1, 
                     (Fy12+Fy13)/m1,
                     k[6],
                     k[7],
                     (Fx21+Fx23)/m2,
                     (Fy21+Fy23)/m2,
                     k[10],
                     k[11],
                     (Fx31+Fx32)/m3,
                     (Fy31+Fy32)/m3])

def rungekutta4(f1, k0, t):
    n = len(t)
    k = np.zeros((n, len(k0)))
    k[0] = k0
    for i in range(n - 1):
        h = t[i+1] - t[i]
        k1 = f1(k[i], t[i]) 
        k2 = f1(k[i] + k1 * h / 2., t[i] + h / 2.)
        k3 = f1(k[i] + k2 * h / 2., t[i] + h / 2.) 
        k4 = f1(k[i] + k3 * h, t[i] + h) 
        k[i+1] = k[i] + (h/6.)*(k1 + 2*k2 + 2*k3 + k4) 
    return k



k0 = np.array([x1, y1, v_x1, v_y1, x2, y2, v_x2, v_y2, x3, y3, v_x3, v_y3])

tf=3600*24*365*70
t4 = np.linspace(0, tf, 50000)
sol4 = rungekutta4(first_order, k0, t4)

plt.plot(sol4[:,0], sol4[:, 1], label='with 21 points')
plt.plot(sol4[:,4], sol4[:, 5], label='with 21 points')
plt.plot(sol4[:,8], sol4[:, 9], label='with 21 points')

#plt.scatter(0,0, s=10)
plt.show()