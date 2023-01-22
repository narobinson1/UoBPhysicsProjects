#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 13:01:40 2022

@author: nrobinson15
"""
import matplotlib.pyplot as plt
import numpy as np
from math import exp

accurate_exp = []
bins = []
bins_index = np.linspace(0,5,num=10)

print(bins_index)
def func1(text_file):
    lines=[]
    f = open(text_file, 'r')
    lines = f.readlines()
    for i in range(0,10):
        y = float(lines[i])
        bins.append(y/800)
        accurate_exp.append(exp(-bins_index[i]))
        


func1("CumulativeMethodData.txt")





plt.plot(bins_index, bins, 'o', ms=4, label='Cumulative Method')
plt.plot(bins_index, accurate_exp, label='$e^{-τ}$')


plt.xlabel("τ")
plt.ylabel("$e^{-τ}$")
plt.legend(prop={'size': 9})


plt.savefig("ExpTauQ1")