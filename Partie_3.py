# -*- coding: utf-8 -*-
"""

Created on Sat May  1 12:55:27 2021

Nom : EL KHALIL

Prénom : Ali

Classe et groupe : 3SC1

Nom du fichier : Partie_3.py

Description : fichier contenant le code pour la résolution de la partie 3 du 
mini projet de Ma322 : équation intégrable.

"""
#-----------------------------------------------
#Importation des modules nécessaires 
#-----------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from math import *

#-----------------------------------------------

h=10**(-3)

#-----------------------------------------------
#Question 2.(a)
#-----------------------------------------------

def Mat(a, b, K, f):
    #initialisation
    N=int((b-a)/h)
    t=np.zeros((N+1,1))
    F=np.zeros((N+1,1))
    A=np.zeros((N+1,N+1))
    #calcul des t_i et les f_i
    for i in range(0,N+1):
        t[i]=a+h*i
        F[i]=f(t[i])
    #calcul des éléments de A
    for i in range(0,N+1):
        A[i,0]=K(t[i],t[0])
        A[i,N]=K(t[i],t[N])
        for j in range(1,N):
            2*K(t[i],t[j])
    #calcul de M
    M=np.eye(N+1)-(h/2)*A
    return t,F,M

def resolution(a,b,K,f):
    t,F,M=Mat(a,b,K,f)
    U=np.linalg.solve(M,F) #Résolution du système MU=F
    return U,t

def comparaison(a,b,K,f,V, figure=False):
    U,t=resolution(a,b,K,f)
    Vtab=V(t)
    if (figure):
        plt.plot(t,U,label="Tracé de U (solution calculée)")
        plt.plot(t,Vtab,label="Tracé de V (solution exacte)")
        plt.xlabel("t")
        plt.ylabel("y(t)")
        plt.title("Courbe des solutions de l'équation intégrable")
        plt.legend()
        plt.show()
    Erreur=np.linalg.norm(U-Vtab)
    return Erreur
    
def K(x,t):
    return 2

def f(x):
    return (np.cos(np.pi*x/2)-8/np.pi)

a=-1
b=1

def V(x):
    return np.cos(np.pi*x/2)


