#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt

i=0

G= 6.67408*(10**-11)
M= 1.989*(10**30)

x0 = 5.2*(10**12)
y0 = 0
v_x0 = 0
v_y0 = 880

k=[] 

def first_order(k, t):
    return np.array([k[2], k[3], -G*M*k[0]/((k[0]**2+k[1]**2)**1.5), -G*M*k[1]/((k[0]**2+k[1]**2)**1.5)])

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

k0 = np.array([x0, y0, v_x0, v_y0])

tf=3600*24*365*76
t4 = np.linspace(0, tf, 50000)
sol4 = rungekutta4(first_order, k0, t4)

plt.plot(sol4[:,0], sol4[:, 1], label='with 21 points')
plt.scatter(0,0, s=10)
plt.show()
