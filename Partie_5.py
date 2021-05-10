# -*- coding: utf-8 -*-
"""
Created on Mon May 10 01:29:21 2021

Nom : EL KHALIL

Prénom : Ali

Classe et groupe : 3SC1

Nom du fichier : Partie_4.py

Description : fichier contenant le code pour la résolution de la partie 5 du 
mini projet de Ma322 : Moteur à Courant Continu (MCC).

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

R=5
L=50*10**-3
Ke=0.2
Kc=0.1
Fm=0.01
Jm=0.05

#condition initiale
Y_0=np.array([0,0])

#-----------------------------------------------
#Question 3.
#-----------------------------------------------

def u(t):
    if t>=10 and t<=50:
        return 5
    else:
        return 0
    
def moteurCC(Y,t):
    i=Y[0]
    w=Y[1]
    Yprim=np.array([(u(t)-R*i-Ke*w)/L,(Kc*i-Fm*w)/Jm])
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
        k1=moteurCC(Y_rk[i-1], t[i-1])
        k2=moteurCC(Y_rk[i-1]+(h/2)*k1,t[i-1]+h/2)
        k3=moteurCC(Y_rk[i-1]+(h/2)*k2,t[i-1]+h/2)
        k4=moteurCC(Y_rk[i-1]+h*k3,t[i-1]+h)
        Y_rk[i]=Y_rk[i-1]+(h/6)*(k1+2*k2+2*k3+k4)
        t.append(h*i)
    if (figure):
        plt.plot(t, Y_rk[:,0])
        plt.title("Tracé de i(t) avec Runge-Kutta 4, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Intensité i (A)")
        plt.show()
        plt.plot(t, Y_rk[:,1])
        plt.title("Tracé de ω(t) avec Runge-Kutta 4, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Vitesse angulaire ω (rad.s^-1)")
        plt.show()
    return (Y_rk,t)


def resolution_odeint(h,tin,tfi,figure=False):
    N=int((tfi-tin)/h)
    t=np.linspace(tin,tfi,N+1)
    Yode =si.odeint(moteurCC, Y_0, t)
    if (figure):
        plt.plot(t, Yode[:,0])
        plt.title("Tracé de i(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Intensité i (A)")
        plt.show()
        plt.plot(t, Yode[:,1])
        plt.title("Tracé de ω(t) avec odeint, pas h = "+str(h))
        plt.xlabel("Temps (sec)")
        plt.ylabel("Vitesse angulaire ω (rad.s^-1)")
        plt.show() 
    return (Yode,t)

def CoupleMoteur(h,tin,tfi,figure=False):
    i,t=(Runge_Kutta_4(h,tin,tfi)[0])[:,0],Runge_Kutta_4(h,tin,tfi)[1]
    Cm=i*Kc
    if (figure):
        plt.plot(t,Cm)
        plt.title("Tracé du couple moteur")
        plt.xlabel("Temps (sec)")
        plt.ylabel("Couple Moteur (N.m)")
        plt.show()

def evolution_vitesse_angulaire(h,tin,tfi,figure=False):
    omega,t=(Runge_Kutta_4(h,tin,tfi)[0])[:,1],Runge_Kutta_4(h,tin,tfi)[1]
    der=[Y_0[1]]
    for i in range (1, len(omega)):
        der.append((omega[i]-omega[i-1])/h)
    if (figure):
        plt.plot(t,der)
        plt.title("Dérivée de la vitesse angulaire")
        plt.xlabel("Temps (sec)")
        plt.ylabel("Acceleration angulaire (rad.s^-2")
        plt.show()
    return der
    
    
h=0.001
tin=0
tfi=80
resolution_odeint(h,tin,tfi,figure=True)
"""Runge_Kutta_4(h,tin,tfi,figure=True)
evolution_vitesse_angulaire(h,tin,tfi,figure=True)
CoupleMoteur(h,tin,tfi,figure=True)"""