#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 17:55:36 2022

@author: nrobinson15
"""

from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



M_Beta025 = np.array([[],[],[],[],[]])
M_Beta033 = np.array([[],[],[],[],[]])
M_Beta066 = np.array([[],[],[],[],[]])
M_Beta074 = np.array([[],[],[],[],[]])

E_Beta025 = np.array([[],[],[],[],[]])
E_Beta033 = np.array([[],[],[],[],[]])
E_Beta066 = np.array([[],[],[],[],[]])
E_Beta074 = np.array([[],[],[],[],[]])



lines=[]
f = open(text_file, 'r')
lines = f.readlines()
     
     



def func1(text_file, beta, n):
    lines=[]
    f = open(text_file, 'r')
    lines = f.readlines()

    for i in range(0,30):
        y = lines[i].split(",")
        M = int(y[0])
        
        E = y[1]
       
        if beta==0.25 and n==1:
            M_Beta025[0] = np.append(M_Beta025[0], M)
            E_Beta025[0] = np.append(E_Beta025[0], M)
       
        if beta==0.25 and n==2:
            M_Beta025[1] = M 
            E_Beta025[1] = E  
        if beta==0.25 and n==3:
            M_Beta025[2] = M
            E_Beta025[2] = E
        if beta==0.25 and n==4:
            M_Beta025[3] = M 
            E_Beta025[3] = E  
        if beta==0.25 and n==5:
            M_Beta025[4] = M  
            E_Beta025[4] = E    
            
        if beta==0.33 and n==1:
            M_Beta033[0] = M   
            E_Beta033[0] = E  
        if beta==0.33 and n==2:
            M_Beta033[1] = M  
            E_Beta033[1] = E  
        if beta==0.33 and n==3:
            M_Beta033[2] = M 
            E_Beta033[2] = E
        if beta==0.33 and n==4:
            M_Beta033[3] = M 
            E_Beta033[3] = E 
        if beta==0.33 and n==5:
            M_Beta033[4] = M  
            E_Beta033[4] = E 
            
        if beta==0.66 and n==1:
            M_Beta066[0] = M
            E_Beta066[0] = E
        if(beta==0.66 and n==2):
            M_Beta066[1] = M  
            E_Beta066[1] = E  
        if(beta==0.66 and n==3):
            M_Beta066[2] = M 
            E_Beta066[2] = E 
        if(beta==0.66 and n==4):
            M_Beta066[3] = M  
            E_Beta066[3] = E  
        if(beta==0.66 and n==5):
            M_Beta066[4] = M
            E_Beta066[4] = E 
            
        if(beta==0.74 and n==1):
            M_Beta074[0] = M  
            E_Beta074[0] = E  
        if(beta==0.74 and n==2):
            M_Beta074[1] = M  
            E_Beta074[1] = E  
        if(beta==0.74 and n==3):
            M_Beta074[2] = M  
            E_Beta074[2] = E 
        if(beta==0.74 and n==4):
            M_Beta074[3] = M   
            E_Beta074[3] = E  
        if(beta==0.74 and n==5):
            M_Beta074[4] = M  
            E_Beta074[4] = E 




              
            
            
func1("1Beta0.25.txt", 0.25, 1)
func1("2Beta0.25.txt", 0.25, 2)
func1("3Beta0.25.txt", 0.25, 3)
func1("4Beta0.25.txt", 0.25, 4)
func1("5Beta0.25.txt", 0.25, 5)

func1("1Beta0.33.txt", 0.33, 1)
func1("2Beta0.33.txt", 0.33, 2)
func1("3Beta0.33.txt", 0.33, 3)
func1("4Beta0.33.txt", 0.33, 4)
func1("5Beta0.33.txt", 0.33, 5)

func1("1Beta0.66.txt", 0.66, 1)
func1("2Beta0.66.txt", 0.66, 2)
func1("3Beta0.66.txt", 0.66, 3)
func1("4Beta0.66.txt", 0.66, 4)
func1("5Beta0.66.txt", 0.66, 5)

func1("1Beta0.74.txt", 0.74, 1)
func1("2Beta0.74.txt", 0.74, 2)
func1("3Beta0.74.txt", 0.74, 3)
func1("4Beta0.74.txt", 0.74, 4)
func1("5Beta0.74.txt", 0.74, 5)

print(M_Beta025)

# # parameters1, covariance1 = curve_fit(func2, R_2, N_2)
# # parameters2, covariance2 = curve_fit(func2, R_3, N_3)
# # parameters3, covariance3 = curve_fit(func2, R_4, N_4)
# # parameters4, covariance4 = curve_fit(func2, R_5, N_5)
# # parameters5, covariance5 = curve_fit(func2, R_6, N_6)

# # err1 = np.sqrt(np.diag(covariance1))
# # err2 = np.sqrt(np.diag(covariance2))
# # err3 = np.sqrt(np.diag(covariance3))
# # err4 = np.sqrt(np.diag(covariance4))
# # err5 = np.sqrt(np.diag(covariance5))


# # array = [parameters1[0],parameters2[0],parameters3[0],parameters4[0],parameters5[0]]

# # print(parameters1)
# # print(parameters2)
# # print(parameters3)
# # print(parameters4)
# # print(parameters5)

# # print(err1)
# # print(err2)
# # print(err3)
# # print(err4)
# # print(err5)


# # print(np.std(array))
# # print(np.mean(array))


# # plt.plot(N_2, R_2, 'o', ms=4,label='s=2, $D_q=1.69 \pm 0.03$')
# # plt.plot(N_3, R_3, 'o', ms=4,label='s=3, $D_q=1.68 \pm 0.05$')
# # plt.plot(N_4, R_4, 'o', ms=4,label='s=4, $D_q=1.62 \pm 0.03$')
# # plt.plot(N_5, R_5, 'o', ms=4,label='s=5, $D_q=2.11 \pm 0.04$')
# # plt.plot(N_6, R_6, 'o', ms=4,label='s=6, $D_q=1.66 \pm 0.04$')
 
 


# # plt.xlabel("ln($\it{r_{max}}$)")
# # plt.ylabel("ln($N_c$)")
# # plt.legend(prop={'size': 9})

# # plt.savefig("FractalDimension1000")