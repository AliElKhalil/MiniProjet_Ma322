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
#On fixe un pas pour toute la partie.
h=10**-3

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
            A[i,j]=2*K(t[i],t[j])
    #calcul de M
    M=np.eye(N+1)-(h/2)*A
    return t,F,M

#-----------------------------------------------
#question 2.(c)
#-----------------------------------------------

def resolution(a,b,K,f, figure=False,titre=False):
    t,F,M=Mat(a,b,K,f)
    U=np.linalg.solve(M,F) #Résolution du système MU=F
    if (figure):
        plt.plot(t,U)
        plt.xlabel("t")
        plt.ylabel("y(t)")
        if (titre):
            plt.title(titre)
        else : 
            plt.title("Courbe de U, solution calculée")
        plt.show()
    return U,t 

#-----------------------------------------------
#question 2.(d)
#-----------------------------------------------

def comparaison(a,b,K,f,V, figure=False, inf=False):
    U,t=resolution(a,b,K,f)
    Vtab=V(t)
    if (figure):
        plt.plot(t,U,label="Tracé de U (solution calculée)", marker='x')
        plt.plot(t,Vtab,label="Tracé de V (solution exacte)")
        plt.xlabel("t")
        plt.ylabel("y(t)")
        plt.title("Courbe des solutions de l'équation intégrable")
        plt.legend()
        plt.show()
    if (inf):
        Erreur=np.linalg.norm(U-Vtab, np.inf)
    else :
        Erreur=np.linalg.norm(U-Vtab)
    return Erreur


if __name__ == '__main__':
    #sous-partie 3.2
    print("Partie 3.2 : Présentation de benchmark utilisée")
    def K(x,t):
        return 2
    def f(x):
        return np.cos(np.pi*x/2)-8/np.pi
    def V(x):
        return np.cos(np.pi*x/2)
    a=-1
    b=1
    resolution(a,b,K,f, figure=True)
    E2=comparaison(a,b,K,f,V,figure=True)
    Einf=comparaison(a,b,K,f,V,inf=True)
    print ("Question 2. (d) :")
    print("L'erreur en norme 2 est : ||U-V||_2 = "+str(E2))
    print("L'erreur en norme infinie est : "+str(Einf))
    print("Fin de la partie 3.2")
    #sous-partie 3.3
    print("Partie 3.3 : Équation de Love en électrostatique")
    def K(x,t):
        return 1/(np.pi*(1+(x-t)**2))
    def f(x):
        return 1
    U,t=resolution(a,b,K,f,figure=True, titre="Courbe de U, solution de l'équation de Love")
    print ("Question 1 :")
    print ("Voici détérminé numériquement le vecteur u : ")
    print(U)
    print("Fin de la partie 3.3")
    print("Fin de la partie 3, merci de votre attention.")

    