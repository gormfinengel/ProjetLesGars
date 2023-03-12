# -*- coding: utf-8 -*-

import sections as sect
#=================== CLASSE MODELE ================================#
#Classe qui reprend tout les instances crée pour modéliser une route et qui va à mettre à jour la résolution numérique du problème

#                ATTRIBUTS
#routes:            tableau avec tout les routes du modèle (ordre sans importance)
#intersection:         tableau avec tout les intersecton du modèle (ordre sans importance)
#carrefours:        tableau avec tout les carrefours du problème (ordre sans importance)
#dt:                 pas de temps du schèma numérique
#h:                    pas d'espace du schèma numérique
#t:                    temps de modélisation
#cellules:            tableau avec toutes les cellules fantome du modèle
#===========================================================================#

class model:
    intersections = None
    routes = None
    dt = None
    h = None
    carrefours = None    
    cellules = None        #Cellules fantomes
    def __init__(self,routes,intersections,carrefours,dt,h, t=0., cellules = None):
         self.intersections = intersections
         self.routes = routes
         self.carrefours = carrefours
         self.cellules = cellules
         self.dt = dt            #pas en temps
         self.h = h                #pas en espace
         self.t = t              #temps en secondes
         
    def resoudre(self):
        #Calcul des flux en entré s'il y a des cellules fantomes
        if self.cellules != None:
            for i in self.cellules:
                i.flux()
            
        
        #Calcul du flux sortie au limites (egale à l'offre s'il y a pas de feux)        
        for i in self.routes:
            if i.sortie:
                if i.feux !=None:
                    if i.feux.passage:
                        i.flux_sortant = sect.demande(i)
                    else:
                        i.flux_sortant = 0
                else:
                    i.flux_sortant = sect.demande(i)
        
        #Calcul flux dans les intersections        
        if self.intersections != None:    
            for i in self.intersections:
                i.flux()
            
        #Calcul flux dans les carrefour
        if self.carrefours != None:
            for i in self.carrefours:
                i.flux()
            
        #Calcul du schema de Godunov
        for i in self.routes:
            sect.godunov(i,self.dt,self.h)
        
        
        #Mise a jour du temps
        self.t=self.t+self.dt
        
