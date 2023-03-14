# -*- coding: utf-8 -*-
import numpy as np
import sections as sect
import intersections as inter
import matplotlib.pyplot as plt

#=================== CLASSE ROUTE_CARREFOUR ================================#
#Classe qui modelise une section de route dans le carrefour
#Classe qui herite de la classe section (voir section.py) avec 3 nouveaux paramètres.

#                ATTRIBUTS
#p_i:         tableau des densites p_i où p_i est la densite des voituires qui veulent sortir a la i-eme sortie
#rentrant:    flux des voitures rentrantes qui vont sortire dans la i-eme sortie 
#sortant:     flux des voitures sortantes qui vont sortire dans la i-eme sortie
#===========================================================================#
class route_carrefour(sect.section):
    
    p_i = None
    rentrant = None
    sortant = None
    
    def __init__(self, nb_sortie, p_max, V_0, T,  nom,  I, U, flux_entrant=0, sortie = False):    
            sect.section.__init__(self, p_max=p_max, V_0=V_0, T=T, nom=nom,  I=I, U=U )
            self.p_i = [[0 for j in range(0,nb_sortie)] for i in self.U]
            self.rentrant = [0 for j in range(0, nb_sortie)]
            self.sortant = [0 for j in range(0,nb_sortie)]
    
    
    #Calcul du nouveau A pour la sortie du carrefour qui dépend des p_i.
    def calcul_A(self, i):
        sors = self.p_i[-1][i]/self.U[-1] if self.U[-1]!=0 else 0
        return [[sors],[1-sors]]
    
    
        

#=================== FONCTION FLUX_I ================================#
#Fonction qui calcule pour une intersection du rond point les flux de l'intersection.

#Parametres:    inter: intersection dans laquelle on va calculer les flux
#        sortie: bolean pour savoir si l'intersection correspondà une sortie du rond point où une entrée
#        A: matrice de passage
#        ne: nombre d'entrée dans le rond point
#        ns: nombre de sortie dans le rond point

#Resultat: modifie les flux rentrant et sortant des routes intermédiaires
#===========================================================================#

def flux_i(inter, sortie,  A ,ne, ns):
    
    #C'est une sortie
    if sortie:
        r1 = inter.route_entree[0]  
        rs = inter.route_sortie[0]
        r2 = inter.route_sortie[1]
    
        #Si il y a pas de densite on arrete
        if int(r1.U[-1]+0.999) == 0:
            r1.sortant = [0 for j in r1.sortant]
            r2.rentrant = [0 for j in r1.sortant]
            return None
        
        f_sortie = rs.flux_entrant
        b_i = [j/r1.U[-1] for j in r1.p_i[-1]]
    
        f_r1 = [j*r1.flux_sortant for j in b_i]
        r1.sortant = f_r1
        
        #Cas ou seulement des voiture qui voulait sortir sortent
        if f_sortie <= f_r1[ns]:
            r2.rentrant = [ r1.sortant[j] - f_sortie*(j==ns)   for j in range(0,len(r1.sortant))]
                    
        #Cas ou des voitures qui voulait sortir dans une autre sortie sortent
        else:       
            f_r1[ns] = 0
            f_sortie = f_sortie - f_r1[ns]
            
            #Distribue la sortie sur les autres voitures
            b_in = [ j/sum(f_r1) for j in f_r1]
            r2.rentrant = [(r1.sortant[j] - b_in[j]*f_sortie)*(ns!=j) for j in range(0,len(r1.sortant))]
    
    #Cas d'une rentree
    else:
            r1 = inter.route_entree[0]  
            re = inter.route_entree[1]
            r2 = inter.route_sortie[0]
    
    #Si il y a pas de densite on arrete
            if r1.U[-1] == 0:
                r1.sortant = [0 for j in r1.sortant]
                r2.rentrant = [re.flux_sortant * A[j][ne] for j in range(0,len(r2.rentrant))]
        
            else:
                b_i = [j/r1.U[-1] for j in r1.p_i[-1]]
                r1.sortant = [j*r1.flux_sortant for j in b_i]
                r2.rentrant = [(re.flux_sortant * A[j][ne]) + r1.sortant[j] for j in range(0,len(r2.rentrant))]
        
     

    
#=================== CLASSE CARREFOUR ================================#
#Classe qui modelise le rond point avec les intersections "classique"
#Classe qui herite de la classe section (voir section.py) avec 3 nouveaux paramètres.

#                ATTRIBUTS
#routes:         Tableau de duplet de sections avec chaque duplet etant une route entrée et une route sortie (represente la double voie) si il y a qu'une voie d'entrée ou qu'une voie que de sortie on fixera l'autre route à None
#            Ordre trigonometrique (sens de la route dans le carrefour)
#A:             Matrice de passage du carrefour entier (taille m*n m étant les entrées et n les sorties)
#taille:        Taille du rond point i.e. le nombre de sortie/entree qu'il y a (se génére automatiquement)
#intersections:     Tableaux des intersections dans le rond point (se remplit automatiquement avec le tableaux routes)
#section_intermediaire:    Tableaux d'instance de la classe route_carrefours qui modélisent les routes interne du rond point
#dt:            Pas de temps du calcul
#h:            Pas en espace
#===========================================================================#    
class carrefour:
    routes = None
    A = None
    taille = None
    intersections = None
    section_intermediaire = None
    dt = None
    h = None
    
    #Constructeur: on passe en paramètres les attributs du rond point et les attributs des routes intermédiaires
    def __init__(self,routes,A,dt,h ,p_max = 150, V_0 = 30., T=.4/3600., longueur = 4,  I= 2):   
            self.routes = routes
            self.A = A
            self.dt = dt
            self.h = h
            self.p_max = p_max
            self.V_0 = V_0
            self.T =T
            self.longueur = longueur
            self.I = I
            #Initialisation valeur taille
            self.taille = len (self.routes)
            taille_sortie = 0
            taille_entree = 0
        
            #Initialisation des tableaux des intersections et des sections de routes dans le carrefour
            self.section_intermediaire = []
            self.intersections = []    
        
            #Calcul du nombre de sortie
            for i in self.routes:
                if i[1]!=None:
                    taille_sortie = taille_sortie + 1
                if i[0]!=None:
                    taille_entree = taille_entree + 1 
            if taille_sortie == 0:
                print("Erreur: pas de sortie")
            if taille_sortie == 0:
                print("Erreur: pas d'entree")
        
            #Initialisation des caracteristiques des sections intermediaires
            N_carr = self.longueur
            u_init_carr = np.zeros(N_carr)
        
            
            #Initialisation des sections intermediaires
            for i in range(0,self.taille):
                self.section_intermediaire.append(route_carrefour(nb_sortie = taille_sortie, p_max=float(self.p_max), V_0=self.V_0, T=self.T,nom="route_carrefour", I=self.I, U = u_init_carr))   
            
            #Initialisation des intersection 
            for i in range(0,self.taille):    
                
                #Si il y a une route de sortie et une route d'entree il faut definir l'intersection interne au carrefour 1-->2: 
                if (self.routes[i][1]!=None and self.routes[i][0]!=None):
                    self.intersections.append(inter.intersection([self.section_intermediaire[i]],[self.routes[i][1],self.section_intermediaire[(i+1)%self.taille]], A=[[0.5],[0.5]],P = [1]))
                                  
                #Si il y a qu'une route de sortie on definit que l'intersection de sortie 1-->2:
                else:
                    if(self.routes[i][1]!=None):
                        self.intersections.append(inter.intersection([self.section_intermediaire[i]],[self.routes[i][1],self.section_intermediaire[(i+1)%self.taille]], A=[[0.5],[0.5]],P = [1]))        
                    
                    #Si il y a qu'une route d'entree on definit que l'intersection d'entree et  2-->1:
                    else:
                        self.intersections.append(inter.intersection([self.section_intermediaire[i],self.routes[i][0]],[self.section_intermediaire[(i+1)%self.taille]] , A= [[1,1]] , P = [1,0]))
    
        
    #Resolution du flux dans tous les intersections du rond point 
    def flux(self):
        nbs = 0
        nbe = 0
        for i in range(0,self.taille):
            re = self.routes[i][0]
            rs = self.routes[i][1]
            r1 = self.section_intermediaire[i]
            r2 = self.section_intermediaire[(i+1)%self.taille]
           
            #Cas ou il y une sortie et un rentree:
            if (re!=None and rs!=None):
            
                #Calcul du A de la sortie:
                A = r1.calcul_A(nbs)
                self.intersections[i].A = A
                
                #Calcul flux avec le nouveau A
                self.intersections[i].flux()
                
                #Calcul des flux_i
                flux_i(self.intersections[i], sortie = True, A= self.A, ne = nbe, ns = nbs)
                
                #Calcul des flux de rentree dans le carrefour s'il y a encore de l'offre
                f_entree = min(sect.demande(re), sect.offre(r2)-r2.flux_entrant)
                re.flux_sortant = f_entree
                r2.flux_entrant = r2.flux_entrant + f_entree
                r2.rentrant = [(re.flux_sortant * self.A[j][nbe]) + r2.rentrant[j] for j in range(0,len(r2.rentrant))]
                
                nbs = nbs+1
                nbe = nbe+1
            
            else:
            #Cas ou il y a qu'une entree:
                if (re!=None):
                    self.intersections[i].flux()
                    #Calcul des flux_i
                    flux_i(self.intersections[i], sortie = False, A = self.A, ne = nbe, ns = nbs)
                    nbe = nbe + 1
                    
                else:
                    #Calcul du A de la sortie:
                    A = r1.calcul_A(nbs)
                    self.intersections[i].A = A
                    #Cas ou il y a qu'une sortie:
                    self.intersections[i].flux()
                    #Calcul des flux_i
                    flux_i(self.intersections[i], sortie = True, A= self.A, ne = nbe, ns = nbs)
                    
                    nbs= nbs + 1
        
        #Godunov sur les routes intermediaires
        for i in range(0,self.taille):
            sect.godunov(self.section_intermediaire[i], self.dt, self.h, carrefour = True)
            
        
    #Affichage des routes intermédiares du rondpoint
    def affichage(self):
        taille = len(self.section_intermediaire[0].p_i[0])
    
        plt.figure(1, figsize= (10,5))       
    
        plt.subplot(131)
        ax = [0,0,0]
        plt.title("route qui relie 1->2")
        for i in range(0,taille):
            l,=plt.plot([j[i] for j in self.section_intermediaire[1].p_i])
            ax[i] = l
        plt.gca().set_xlabel('x(m)')
        plt.gca().xaxis.set_ticklabels([ '0', '10', '20','30', '40', '50','60', '70','80', '90', '100'], rotation = 0,  fontsize = 10)
        plt.gca().set_ylabel('densite')
        plt.legend([ax[0], ax[1], ax[2]],["sortie 1", "sortie 2", "sortie 3"],loc=1, prop={'size': 10}) 
    
        plt.subplot(132) 
        plt.title("route qui relie 2->3")
        for i in range(0,taille):
            plt.plot([j[i] for j in self.section_intermediaire[2].p_i])
        plt.gca().xaxis.set_ticklabels([ '0', '10', '20','30', '40', '50','60', '70','80', '90', '100'], rotation = 0,  fontsize = 10)
        plt.gca().set_xlabel('x(m)')
            
        plt.subplot(133)    
        plt.title("route qui relie 3->1")
        for i in range(0,taille):
            plt.plot([j[i] for j in self.section_intermediaire[0].p_i])
        plt.gca().xaxis.set_ticklabels([ '0', '10', '20','30', '40', '50','60', '70','80', '90', '100'], rotation = 0,  fontsize = 10)
        plt.gca().set_xlabel('x(m)')
        
        plt.suptitle("T = 00mins 30 sec")
        plt.show()
        
    
    
