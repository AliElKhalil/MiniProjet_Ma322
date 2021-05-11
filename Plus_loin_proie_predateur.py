# -*- coding: utf-8 -*-
"""
Created on Mon May 10 12:56:01 2021

Nom : EL KHALIL

Prénom : Ali

Classe et groupe : 3SC1

Nom du fichier : Partie_7.py

Description : fichier contenant le code pour la résolution de la partie 6 du 
mini projet de Ma322 : Modèle proie-prédateur.

"""

#-----------------------------------------------
#Importation des modules nécessaires 
#-----------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from math import *
import scipy.integrate as si

#-----------------------------------------------
#Condition initiale
Y_0=np.array([5,3,4])
alpha1=3
beta1=1
alpha2=2
beta2=3
alpha3=3
beta3=1
beta4=1
Tin=0
Tfi=20
h=10**-4

#-----------------------------------------------
#Question 3.
#-----------------------------------------------

def ProiePredateur(Y, t) :
    Yprime = np.zeros(3)
    Yprime[0] = alpha1*Y[0]-beta1*Y[0]*Y[1]
    Yprime[1] = -alpha2*Y[1]+beta2*Y[0]*Y[1]-beta3*Y[1]*Y[2]
    Yprime[2] = -alpha3*Y[2]+beta4*Y[1]*Y[2]
    return Yprime

def Euler(h, Tin, Tfi, figure=False):
    N=int ((Tfi-Tin)/h)
    Y_e=np.zeros((N+1,3))
    Y_e[0]=Y_0
    t=[Tin]
    for i in range(1,N+1):
        Y_e[i]=Y_e[i-1]+h*ProiePredateur(Y_e[i-1], t[i-1])
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_e[:,0],label="Proie")
        plt.plot(t, Y_e[:,1],label="Prédateur 1")
        plt.plot(t, Y_e[:,2],label="Prédateur 2")
        plt.title("Modèle proie-prédateur Euler")
        plt.xlabel("Temps (u.t.)")
        plt.ylabel("Nombre d'animal")
        plt.legend()
        plt.show()
    return (Y_e)

#-----------------------------------------------
#Question 4.
#-----------------------------------------------

def resolution_odeint(h,tin,tfi,figure=False):
    N=int ((tfi-tin)/h)
    N=int((tfi-tin)/h)
    t=np.linspace(tin,tfi,N+1)
    Yode =si.odeint(ProiePredateur, Y_0, t)
    if (figure):
        plt.plot(t, Yode[:,0],label="Proie")
        plt.plot(t, Yode[:,1],label="Prédateur 1")
        plt.plot(t, Yode[:,2],label="Prédateur 2")
        plt.title("Modèle proie-prédateur Runge-Kutta 4")
        plt.xlabel("Temps (u.t.)")
        plt.ylabel("Nombre d'animal")
        plt.legend()
        plt.show()
    return (Yode)
    
resolution_odeint(h,Tin, Tfi,figure=True)  