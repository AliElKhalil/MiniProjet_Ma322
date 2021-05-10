# -*- coding: utf-8 -*-
"""

Created on Sat May  4 10:59:30 2021

Nom : EL KHALIL

Prénom : Ali

Classe et groupe : 3SC1

Nom du fichier : Partie_4.py

Description : fichier contenant le code pour la résolution de la partie 4 du 
mini projet de Ma322 : Circuit RLC.

"""
#-----------------------------------------------
#Importation des modules nécessaires 
#-----------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from math import *
import scipy.integrate as si

#-----------------------------------------------
#on fixe les constantes
C=10**-6
R=3
L=0.5

t0=0
tf=2
#condition initiale
Y_0=np.array([0,0])
#-----------------------------------------------
#Question 3.
#-----------------------------------------------

def e(t):
    if (t>=0):
        return 10
    else:
        return 0


def rlcprim(Y,t):
    s=Y[0]
    i=Y[1]
    Yprim=np.array([i/C,(e(t)-s-R*i)/L])
    return Yprim

#-----------------------------------------------
#Question 4 et 5
#-----------------------------------------------

def Runge_Kutta_4(h,tin,tfi,figure=False):
    N=int ((tfi-tin)/h)
    Y_rk=np.zeros((N+1,2))
    Y_rk[0]=Y_0
    t=[tin]
    for i in range(1, N+1):
        k1=rlcprim(Y_rk[i-1], t[i-1])
        k2=rlcprim(Y_rk[i-1]+(h/2)*k1,t[i-1]+h/2)
        k3=rlcprim(Y_rk[i-1]+(h/2)*k2,t[i-1]+h/2)
        k4=rlcprim(Y_rk[i-1]+h*k3,t[i-1]+h)
        Y_rk[i]=Y_rk[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_rk[:,0])
        plt.title("Tracé de s(t) avec Runge-Kutta 4, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Tension s (V)")
        plt.show()
        plt.plot(t, Y_rk[:,1])
        plt.title("Tracé de i(t) avec Runge-Kutta 4, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Intensité i (A)")
        plt.show()
    return (Y_rk,t)


def resolution_odeint(h,tin,tfi,figure=False):
    N=int((tfi-tin)/h)
    t=np.linspace(tin,tfi,N+1)
    Yode =si.odeint(rlcprim, Y_0, t)
    if (figure):
        plt.plot(t, Yode[:,0])
        plt.title("Tracé de s(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Tension s (V)")
        plt.show()
        plt.plot(t, Yode[:,1])
        plt.title("Tracé de i(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Intensité i (A)")
        plt.show() 
    return (Yode,t)
    

if __name__ == '__main__':
    c=1
    while(c==1): 
        print("Partie 4 : circuit RLC.")
        h=float0.001(input("Saisir le pas souhaité : "))
        Runge_Kutta_4(h,t0,tf,figure=True)
        resolution_odeint(h,t0,tf,figure=True)
        c=2
        while(c!=1 and c!=0):
            c=int(input("Voulez vous recommencer ? 1 si oui, 0 si non : "))
    print("Fin de la partie 4., merci pour votre attention.")
    
    