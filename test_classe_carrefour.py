# -*- coding: utf-8 -*-
import sys
sys.path.append("../Model/LWR")

import numpy as np
import sections as sect

import carrefour as carr
import modele as mod
import matplotlib.pyplot as plt


#Initialisation du modèle route
N=10                     #Longueur du modèle
u_init1 = np.zeros(N)
u_init2 = np.zeros(N)
u_init3 = np.zeros(N)
u_init4 = np.zeros(N)
u_init1[10:50] = 0
u_init2[30:80] = 0
u_init3[20:90] = 0
h=0.02                  #taille de chaque maille (km)
vmax= 50.               #Vitesse maximale du reseau routier (km/h)
dt = 0.8*(h/vmax)      #Pas de tems (h)

#Initialisation du modèle du carrefour
V_0 = vmax              #Vitesse souhaite (km/h)



#Definition routes entrantes sortantes
route1 = sect.section(p_max=100., V_0=V_0, T=1.4/3600,nom="route1",I = 1., U = u_init1, flux_entrant = 2000) 
route2 = sect.section(p_max=100., V_0=V_0, T=1.4/3600,nom="route2",I = 1., U = u_init2, flux_entrant = 1500) 
route3 = sect.section(p_max=100., V_0=V_0, T=1.4/3600,nom="route3",I = 1., U = u_init3, flux_entrant = 1800) 
route4 = sect.section(p_max=100., V_0=V_0, T=1.4/3600,nom="route4",I = 1., U = u_init4, sortie = True) 
route5 = sect.section(p_max=100., V_0=V_0, T=1.4/3600,nom="route5",I = 1., U = u_init4, sortie = True) 
route6 = sect.section(p_max=100., V_0=V_0, T=1.4/3600,nom="route6",I = 1., U = u_init4, sortie = True)


#A= [[0,0,1],[0,0,0],[1,1,0]]
A= [[1./3,1./3,1./3],[1./3,1./3,1./3],[1./3,1./3,1./3]]
carrefour1=carr.carrefour([[route1,route4],[route2,route5],[route3,route6]], A= A, dt = dt, h =h,p_max = 100., V_0 = 50., T=1./3600., longueur = 4,  I= 1)

#Modele
modele = mod.model(routes = [route1,route2,route3,route4,route5,route6],intersections = None, carrefours = [carrefour1], dt =dt, h=h)

iteration = 0
Sauve1 = []
Sauve2 = []
Sauve3 = []

while iteration<4000:
    
    #résout le modèle
    modele.resoudre()

    #Calcul du temps en heures:minutes:secondes
    t=modele.t*3600
    heure = int(t/3600)
    mins = int((t-heure*3600)/60)
    sec = int(t-heure* 3600 - mins * 60)
    iteration = iteration +1
    #Affichage tous les 15 secondes
    
    if (sec+mins*60)%15 == 0:
        print ("########################################################   FIGURE    ######################################################################")
        #Affiche les routes intermédiaires
        carrefour1.affichage()
        
        #Affiche les routes d'entree/sortie
        plt.figure(1, figsize= (14,10))       
        plt.subplot(231)
        plt.title("route d'entree 1")
        plt.gca().set_ylabel('Densite')
        plt.gca().set_ylim(-1,200)
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '200', '400', '600', '800','1000'], rotation = 0,  fontsize = 10)
        plt.cla    
        plt.plot(route1.U)
           
        plt.subplot(232)
        plt.title("route d'entree 2")
        plt.gca().set_ylim(-1,100)
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '200', '400', '600', '800','1000'], rotation = 0,  fontsize = 10)
        plt.cla    
        plt.plot(route2.U)
    
        plt.subplot(233)
        plt.title("route d'entree 3")
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '200', '400', '600', '800','1000'], rotation = 0,  fontsize = 10)
        plt.gca().set_ylim(-1,200)          
        plt.cla    
        plt.plot(route3.U)
              
        plt.subplot(234)
        plt.title("route de sortie 1")
        plt.gca().set_ylim(-1,100)    
        plt.gca().set_ylabel('Densite')            
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '200', '400', '600', '800','1000'], rotation = 0,  fontsize = 10)
        plt.cla    
        plt.plot(route4.U)
    
    
        plt.subplot(235)
        plt.title("route de sortie 2")
        plt.gca().set_ylim(-1,100)          
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '200', '400', '600', '800','1000'], rotation = 0,  fontsize = 10)
        plt.cla    
        plt.plot(route5.U)
    
        plt.subplot(236)
        plt.title("route de sortie 3")
        plt.gca().set_ylim(-1,100)    
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '200', '400', '600', '800','1000'], rotation = 0,  fontsize = 10)
        plt.cla    
        plt.plot(route6.U)
    
        plt.suptitle("T = " +str(mins) + "mins " + str(sec) +"sec")
        plt.show()
        
        


