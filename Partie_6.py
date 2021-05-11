# -*- coding: utf-8 -*-
"""
Created on Mon May 10 11:13:27 2021

Nom : EL KHALIL

Prénom : Ali

Classe et groupe : 3SC1

Nom du fichier : Partie_6.py

Description : fichier contenant le code pour la résolution de la partie 6 du 
mini projet de Ma322 : Mouvement d'une fusée.

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
a0=8*10**3
g=9.81
k0=0.1
u=2*10**3

tin=0
tfi=100

#condition initiale
Y_0=np.array([0,400,0])

#-----------------------------------------------
#Question 1.
#-----------------------------------------------

def fusee(Y, t) :
    D = 4
    Yprime = np.zeros(3)
    if (Y [1] < 80) :
        Y [1] = 80
        D = 0
    Yprime[0] = D*u/Y[1]-g-k0*np.exp(-Y[2]/a0)*(Y[0]**2)/Y[1]
    Yprime[1] = -D
    Yprime[2] = Y[0]
    return Yprime

#-----------------------------------------------
#Question 2.
#-----------------------------------------------

def resolution_odeint(h,tin,tfi,figure=False):
    N=int((tfi-tin)/h)
    t=np.linspace(tin,tfi,N+1)
    Yode =si.odeint(fusee, Y_0, t)
    if (figure):
        plt.plot(t, Yode[:,0])
        plt.title("Tracé de v(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Vitesse v (m.s^-1)")
        plt.show()
        plt.plot(t, Yode[:,1],label="m(t)")
        plt.plot(np.linspace(tin,tfi,25),80*np.ones((25)),label="y = 80")
        plt.title("Tracé de m(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("masse m (kg)")
        plt.legend()
        plt.show()
        plt.plot(t, Yode[:,2])
        plt.title("Tracé de z(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("altitude z (m)")
        plt.show()
    return (Yode,t)


if __name__ == '__main__':
    c=1
    while(c==1): 
        print("Partie 6 : Lancement d'une fusée.")
        h=float(input("Saisir le pas souhaité : "))
        resolution_odeint(h,tin,tfi,figure=True)
        c=2
        while(c!=1 and c!=0):
            c=int(input("Voulez vous recommencer ? 1 si oui, 0 si non : "))
    print("Fin de la partie 4., merci pour votre attention.")