import sys

import modele as mod
import sections as secc

# -*- coding: utf-8 -*-
import sys
sys.path.append("Model/LWR")

import numpy as np
import sections as sect
import intersections as inter
import modele as mod
import matplotlib.pyplot as plt


#Initialisation du modèle
N=20                            #Longueur du modèle (en maille)
u_init1 = np.zeros(N)               
u_init2 = np.zeros(N)
u_init1[int(N/2):int(3*N/4)] = 100        
h=0.1                           #taille de chaque maille (km)
vmax= 50.                       #Vitesse maximale du systeme routier
dt = 0.99*(h/vmax)              #Condition CFL
print ("dt = " + str(dt))  

#Definition routes
route1 = sect.section(p_max=100., V_0=50., T=1.4/3600,nom="1", I = 1., U = u_init1) 
route2 = sect.section(p_max=100., V_0=50., T=1.4/3600,nom="2", I = 1., U = u_init1) 
route3 = sect.section(p_max=100., V_0=50., T=1.4/3600,nom="3", I = 1., U = u_init1) 
route4 = sect.section(p_max=100., V_0=50., T=1.4/3600,nom="4", I = 1., U = u_init2 , sortie = True) 
route5 = sect.section(p_max=100., V_0=50., T=1.4/3600,nom="5", I = 1., U = u_init2, sortie = True) 
route6 = sect.section(p_max=100., V_0=50., T=1.4/3600,nom="6", I = 1., U = u_init2, sortie = True) 


#Definition intersection
intersection1 = inter.intersection([route1,route2],[route3,route4],A=[[0.5,0.5],[0.5,0.5]],P = [0.3,0.7])
intersection2 = inter.intersection([route3],[route5,route6],A=[[0.4],[0.6]],P = [1])
intersection3 = inter.intersection([route5],[route1,route2],A=[[0.2],[0.8]],P = [1])


#Modele
modele = mod.model([route1,route2,route3,route4,route5,route6],[intersection1,intersection2,intersection3],None,dt,h)


iteration = 1

while iteration<100:
    #Resoudre le modele
    modele.resoudre()
    
    iteration+=1
    
    #Plot tous les 5 iterations
    if iteration%5 ==0:
        print ("##############   FIGURE    ##############")
        plt.figure()       
        plt.subplot(331)
        plt.gca().set_ylim(-1,100)
        plt.cla    
        plt.plot(route1.U)
           
        plt.subplot(332)
        plt.gca().set_ylim(-1,100)
        plt.cla    
        plt.plot(route2.U)
              
        plt.subplot(333)
        plt.gca().set_ylim(-1,100)
        plt.cla    
        plt.plot(route3.U)
              
        plt.subplot(334)
        plt.gca().set_ylim(-1,100)
        plt.cla    
        plt.plot(route4.U)
        
        
        plt.subplot(335)
        plt.gca().set_ylim(-1,100)
        plt.cla    
        plt.plot(route5.U)
        
        plt.subplot(336)
        plt.gca().set_ylim(-1,100)
        plt.cla    
        plt.plot(route6.U)
        plt.show()

