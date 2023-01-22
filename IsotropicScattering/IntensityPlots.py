#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 13:01:40 2022

@author: nrobinson15
"""
import matplotlib.pyplot as plt
import numpy as np
from math import exp

Isotropic_intensity = []

bins = []
bins_index = np.linspace(0,9,num=10)


def func1(text_file):
    lines=[]
    f = open(text_file, 'r')
    lines = f.readlines()
    for i in range(0,10):
        y = float(lines[i])
        Isotropic_intensity.append(y)
       
            
        


func1("IsotropicScattering.txt")





plt.plot(bins_index, Isotropic_intensity, 'o', ms=4, label='Isotropic Scattering')



plt.xlabel("Bin number")
plt.ylabel("Relative intensity")
plt.legend(prop={'size': 9})


plt.savefig("ScatteringExperiment")