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
Y_0=np.array([5,3])
alpha1=3
beta1=1
alpha2=2
beta2=1
Tin=0
Tfi=20
h=10**-4
#-----------------------------------------------
#Question 1.
#-----------------------------------------------

def y1_SansPredateur(t):
    return Y_0[0]*np.exp(alpha1*t)


def trace_y1_SansPredateur(Tin,Tfi):
    t=np.linspace(Tin,Tfi,1000)
    Y1=y1_SansPredateur(t)
    plt.plot(t,Y1)
    plt.title("Courbe d'évolution des proies en l'absence de prédateur")
    plt.xlabel("temps t (u.t.)")
    plt.ylabel("Nombre d'animal")
    plt.show()
    
#-----------------------------------------------
#Question 2.
#-----------------------------------------------

def y2_SansProie(t):
    return Y_0[1]*np.exp(-alpha2*t)


def trace_y2_SansProie(Tin,Tfi):
    t=np.linspace(Tin,Tfi,1000)
    Y2=y2_SansProie(t)
    plt.plot(t,Y2)
    plt.title("Courbe d'évolution des prédateur en l'absence de proie")
    plt.xlabel("temps t (u.t.)")
    plt.ylabel("Nombre d'animal")
    plt.show()
    
#-----------------------------------------------
#Question 3.
#-----------------------------------------------

def ProiePredateur(Y, t) :
    Yprime = np.zeros(2)
    Yprime[0] = alpha1*Y[0]-beta1*Y[0]*Y[1]
    Yprime[1] = -alpha2*Y[1]+beta2*Y[0]*Y[1]
    return Yprime

def Euler(h, Tin, Tfi, figure=False):
    N=int ((Tfi-Tin)/h)
    Y_e=np.zeros((N+1,2))
    Y_e[0]=Y_0
    t=[Tin]
    for i in range(1,N+1):
        Y_e[i]=Y_e[i-1]+h*ProiePredateur(Y_e[i-1], t[i-1])
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_e[:,0],label="Proie")
        plt.plot(t, Y_e[:,1],label="Prédateur")
        plt.title("Modèle proie-prédateur Euler")
        plt.xlabel("Temps (u.t.)")
        plt.ylabel("Nombre d'animal")
        plt.legend()
        plt.show()
    return (Y_e)

#-----------------------------------------------
#Question 4.
#-----------------------------------------------

def Runge_Kutta_4(h,tin,tfi,figure=False):
    N=int ((tfi-tin)/h)
    Y_rk=np.zeros((N+1,2))
    Y_rk[0]=Y_0
    t=[tin]
    for i in range(1, N+1):
        k1=ProiePredateur(Y_rk[i-1], t[i-1])
        k2=ProiePredateur(Y_rk[i-1]+(h/2)*k1,t[i-1]+h/2)
        k3=ProiePredateur(Y_rk[i-1]+(h/2)*k2,t[i-1]+h/2)
        k4=ProiePredateur(Y_rk[i-1]+h*k3,t[i-1]+h)
        Y_rk[i]=Y_rk[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_rk[:,0],label="Proie")
        plt.plot(t, Y_rk[:,1],label="Prédateur")
        plt.title("Modèle proie-prédateur Runge-Kutta 4")
        plt.xlabel("Temps (u.t.)")
        plt.ylabel("Nombre d'animal")
        plt.legend()
        plt.show()
    return (Y_rk)

#-----------------------------------------------
#Question 5.
#-----------------------------------------------

def portrait_de_phase(h,Tin,Tfi):
    Y_e=Euler(h, Tin, Tfi)
    Y_rk=Runge_Kutta_4(h, Tin, Tfi)
    plt.plot(Y_e[:,0],Y_e[:,1])
    plt.title("Portrait de phase avec la méthode d'Euler")
    plt.xlabel("y1(t)")
    plt.ylabel("y2(t)")
    plt.show()
    plt.plot(Y_rk[:,0],Y_rk[:,1])
    plt.title("Portrait de phase avec la méthode de Runge-Kutta d'ordre 4")
    plt.xlabel("y1(t)")
    plt.ylabel("y2(t)")
    plt.show()
  
    

if __name__ == '__main__':
    c=1
    while(c==1): 
        print("Partie 7 : Modèle proie prédateur.")
        h=float(input("Saisir le pas souhaité : "))
        alpha1=float(input("saisir le paramètre alpha1 : "))
        alpha2=float(input("saisir le paramètre alpha2 : "))
        beta1=float(input("Saisir le paramètre beta1 : "))
        beta2=float(input("Saisir le paramètre beta2 : "))
        Y_0[0]=int(input("Saisir le nombre initial de proie : "))
        Y_0[1]=int(input("Saisir le nombre initial de prédateur : "))
        Euler(h,Tin,Tfi,figure=True)
        Runge_Kutta_4(h, Tin, Tfi,figure=True)
        portrait_de_phase(h, Tin, Tfi)
        c=2
        while(c!=1 and c!=0):
            c=int(input("Voulez vous recommencer ? 1 si oui, 0 si non : "))
    print("Fin de la partie 7., merci pour votre attention.")