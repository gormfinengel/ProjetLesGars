#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 14:28:58 2018

@author: abernabe
"""

import matplotlib.pyplot as plt
import modele as mod
import numpy as np
import sections as s

#Caractéristique de la route
N = 250                 #Taille de la route (en maille)
U_tot = np.zeros(N)     #Initialisation de la route
U_tot[0:30] = 20
p_c = 100.              #Densite critique (veh/km)
p_max = 150.            #Densite maximale (veh/km)

#Caractéristique conducteurs
V_0 = 110.            #Vitesse souhaitée par chaque groupe de conducteurs
V_c = 120.                             #Vitesse critique (pour l'offre de la route)


h = 1.                      #taille de chaque maille (km)
dt =  0.99*(h/V_0)     #CFL

#Definition de la route multi classe
route = s.section(p_max=p_max, V_0=V_0, p_c = p_c, nom="Route test", I = 1, U = U_tot)

#Definition du modele multi classe
model_r = mod.model([route], intersections= None, carrefours= None, dt= dt,h= h)


iteration = 0

Sauve = [] #Tableau pour sauvegarder

while iteration < 200:
    model_r.resoudre()
    
    #Plot
    plt.plot(route.U, color = 'green')
    
    #Sauvegarde l'état de la route
    Sauve.append(np.copy(route.U))
    iteration = iteration +1
    
    

#Plot
fig = plt.figure(figsize= (10,5))
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.imshow(Sauve, interpolation='nearest', cmap=plt.cm.ocean, aspect='auto' ,vmin=0,vmax=30)
plt.colorbar()
plt.show()
