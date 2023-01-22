#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 18:10:38 2022

@author: nrobinson15
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
i=0

G= 6.67408*(10**-11)
M= 1.989*(10**30)

x1_0 = 2.52*1.496*10**11
y1_0 = 0
v_x1_0 = 0
v_y1_0 = sqrt(G*M/x1_0)

x2_0 = 5.24*1.496*10**11
y2_0 = 0
v_x2_0 = 0
v_y2_0 = sqrt(G*M/x2_0)


k=[] 


def first_order(k, t):
    return np.array([k[2], 
                     k[3], 
                     -G*M*k[0]/((k[0]**2+k[1]**2)**1.5), 
                     -G*M*k[1]/((k[0]**2+k[1]**2)**1.5), 
                     k[6], 
                     k[7], 
                     -G*M*k[4]/((k[4]**2+k[5]**2)**1.5), 
                     -G*M*k[5]/((k[4]**2+k[5]**2)**1.5)])

def rungekutta4(f, k0, t):
    n = len(t)
    k = np.zeros((n, len(k0)))
    k[0] = k0
    for i in range(n - 1):
        
        h = t[i+1] - t[i]
        k1 = f(k[i], t[i]) 
        k2 = f(k[i] + k1 * h / 2., t[i] + h / 2.)
        k3 = f(k[i] + k2 * h / 2., t[i] + h / 2.) 
        k4 = f(k[i] + k3 * h, t[i] + h) 
        k[i+1] = k[i] + (h/6.)*(k1 + 2*k2 + 2*k3 + k4) 
        
    return k

print(sqrt(G*(10**-3)*M/x1_0))
print(sqrt(G*4*(10**-2)*M/x2_0))

k0 = np.array([x1_0, y1_0, v_x1_0, v_y1_0, x2_0, y2_0, v_x2_0, v_y2_0])
print(k0)
tf=3600*24*365*29
t4 = np.linspace(0, tf, 50000)
sol4 = rungekutta4(first_order, k0, t4)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')

plt.plot(sol4[:,0], sol4[:, 1], label='with 21 points')
plt.plot(sol4[:,4], sol4[:,5])


plt.scatter(0,0, s=10)
plt.show()



