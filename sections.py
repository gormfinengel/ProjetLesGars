# -*- coding: utf-8 -*-
import numpy as np
from functools import partial



#=================  GESTION DES MODELES  ================================#
#case == 0 => modele triangulaire
#case == 1 => modele greenshields
#case == 2 => modele greenshields generalisé
case = 0
n_green = 4.
#===========================================================================#

#=================  FONCTION RIEMANN  ================================#
#Résout le problème de Riemann élémentaire. Fonctionne lorsque la fonction flux est concave

#paramètres:    mod: section de route ou l'on calcule la solution de Rieman
#        p_L: densité à gauche
#        p_R: densité à droite
#===========================================================================#

def Riemann(mod,p_L,p_R):
    if (p_L < p_R):
        if (mod.f(p_L) < mod.f(p_R)):
            return p_L
        else:
            return p_R
    elif (p_L > p_R):
        if (0 < mod.df(p_L)):
            return p_L
        elif (0 > mod.df(p_R)):
            return p_R
        else:
            return mod.finv(0)
    else:
        return p_L
Riemann = np.vectorize(Riemann)
#=================  FONCTION OFFRE/DEMANDE  ================================#
#Calcule l'offre/demande d'une route

#paramètres:    mod: section de route ou l'on calcule l'offre/demande

#Remarque: on prend en compte les feux pour la demande (un feux rouge met automatiquement la demande à 0)
#===========================================================================#
def offre(mod):
    if (mod.U[0] < mod.p_c):
        return mod.f(mod.p_c)
    else: 
        return mod.f(mod.U[0])
        
def demande(mod):
    if mod.feux == None:
        vert = True
    else:
        vert = mod.feux.passage()
    if vert:
        if (mod.U[-1] < mod.p_c):
            resultat = mod.f(mod.U[-1]) #- mod.correction
            if resultat < 0:
                return 0
            return resultat
        else:
            return mod.f(mod.p_c) #- mod.correction
    else:
        return 0

#=================  FONCTION GODUNOV  ================================#
#Résolutiion du schèma numérique de Godunov sur la densité.
#Il calcule aussi  le schéma numérique sur les densité partiel p_i pour les sections de route intermédiaires des rond point.
#Met à jours le temps d'attente moyen

#paramètres:    mod: section de route  mise à jour
#        dt: pas en temps du schèma numérique
#        h:pas en espace du schèma numérique
#        carrefour: boolean pour savoir si la mod est une route intermédiaire dans un rond point (pour mettre à jours les densité partiels p_i)
#===========================================================================#
def godunov(mod, dt, h, carrefour = False):
    
    N = int(len(mod.U))


    #Mise à jour de la densité dans chaque maille
    Un = np.copy(mod.U)
    mod.U[1:N-1] = Un[1:N-1] - (dt/h)*(mod.f(Riemann(mod,Un[1:N-1],np.roll(Un,-1)[1:N-1])) - mod.f(Riemann(mod,np.roll(Un,1)[1:N-1],Un[1:N-1])))
    mod.U[0]=Un[0]- (dt/h)*(mod.f(Riemann(mod,Un[0],Un[1]))-mod.flux_entrant)
    mod.U[N-1]=Un[N-1] - (dt/h)*(mod.flux_sortant - mod.f(Riemann(mod,Un[N-2],Un[N-1])))    
    
#Dev: Pratik (begin)    
#    Un = np.copy(mod.U)
#    mod.U[[1:N-1],[1:N-1]] = Un[[1:N-1],[1:N-1]] - (dt/h)*(mod.f(Riemann(mod,Un[1:N-1],np.roll(Un,-1)[1:N-1])) - mod.f(Riemann(mod,np.roll(Un,1)[1:N-1],Un[1:N-1])))
#    mod.U[0]=Un[0]- (dt/h)*(mod.f(Riemann(mod,Un[0],Un[1]))-mod.flux_entrant)
#    mod.U[N-1]=Un[N-1] - (dt/h)*(mod.flux_sortant - mod.f(Riemann(mod,Un[N-2],Un[N-1])))    
    
    #Mise à jour du nombre cumulatif de vehicule
    n = np.copy(mod.N)
    mod.N[1:N] = n[1:N] + dt * mod.f(Riemann(mod,np.roll(Un,1)[1:N],Un[1:N]))
    mod.N[0] = n[0] + dt * mod.flux_entrant
    
    #Mise à jour du temps d'attente 
    nb_v = [i*h for i in Un]
    t_old = np.copy(mod.t_att)
    mod.t_att[1:N-1] = t_old[1:N-1]*nb_v[1:N-1] - dt*(mod.f(Riemann(mod,Un[1:N-1],np.roll(Un,-1)[1:N-1]))*t_old[1:N-1]  \
                                                      - mod.f(Riemann(mod,np.roll(Un,1)[1:N-1],Un[1:N-1]))*np.roll(t_old,1)[1:N-1])
    mod.t_att[0] = t_old[0]*nb_v[0]  - dt*(mod.f(Riemann(mod,Un[0],Un[1]))*t_old[0] - mod.flux_t_entrant)
    mod.t_att[N-1] = t_old[N-1]*nb_v[N-1] - dt*(mod.flux_sortant * t_old[N-1] - mod.f(Riemann(mod,Un[N-2],Un[N-1]))* t_old[N-2])
    mod.t_att =[mod.t_att[i]/(h*mod.U[i])+dt if h*mod.U[i]>0.01 else 0 for i in range(len(nb_v))] 

    #Mise a jour du temps du feux
    if mod.feux != None:
        mod.feux.t = mod.feux.t + mod.feux.dt
    
    #Calcul des p_i pour les sections du carrefour
    if carrefour:
        old_pi = np.copy(mod.p_i)    
        for i in range(0,len(mod.p_i)):
            for j in range(0,len(mod.p_i[0])):        
            #Calcul premier element
                if i == 0:
                    b_i = old_pi[i][j]/Un[i] if Un[i]!=0 else 0    
                    mod.p_i[i][j] = old_pi[i][j] - (dt/h)*(b_i*mod.f(Riemann(mod,Un[0],Un[1])) - mod.rentrant[j])                
                #Calcul dernier element
                elif i == len(mod.p_i)-1:
                    b_i1 = old_pi[i-1][j]/Un[i-1] if Un[i]!=0 else 0
                    mod.p_i[i][j] = old_pi[i][j] - (dt/h)*(mod.sortant[j] - b_i1 * mod.f(Riemann(mod,Un[i-1],Un[i])))        
                #Calcul du reste    
                else:
                    b_i1 = old_pi[i-1][j]/Un[i-1] if Un[i]!=0 else 0
                    b_i = old_pi[i][j]/Un[i] if Un[i]!= 0 else 0
                    mod.p_i[i][j] = old_pi[i][j] - (dt/h)*(b_i * mod.f(Riemann(mod,Un[i],Un[i+1])) - b_i1 * mod.f(Riemann(mod,Un[i-1],Un[i])))
            
    
    
#================  MODELE TRIANGULAIRE  ===============================#
if case == 0:
    #Fonction vitesse    #vitesse triangulaire
    def V(p_c,T,V_0,p_max,l_eff,p):
        if (p <= p_c):
            return V_0
        elif (p_c < p <= p_max ):
            return 1/(T * p)-l_eff/T
        else:
            print("ein problem")
            print("densite ="  + str(p))
            print("densite_max" + str(p_max))
            return 
    V = np.vectorize(V)
    
    #Fonction flux TRIANGULAIRE
    def F(p_c,T,l_eff,V_0,p):
        if (p <= p_c):
            return V_0*p
        else:
            return (1./T)*(1-p*l_eff)
    
    F = np.vectorize(F)
    
    #Fonction dérivée du flux
    def DF(p_c,V_0,T,l_eff,p):    
        if (p < p_c):
            return V_0
        else:
            return -l_eff/T  
    DF = np.vectorize(DF)
    
    #Fonction inverse de la dérivée du flux
    def FINV(p_c,p_max,p):
        return p_c+p_max*0
    FINV = np.vectorize(FINV)        
#===================  MODELE GREENSHIELDS  ======================================#
if case == 1:
    n_green = 1
    def V(p_c,T,V_0,p_max,l_eff,p):   #VITESSE GREENSHIELDS
        return V_0*(1 - p/p_max)
    V = np.vectorize(V)
    
        
    def F(p_c,T,l_eff,V_0,p):       #FLUX GRRENSHIELDS
        return p*V_0*(1 - p*l_eff)
    F = np.vectorize(F)
    
    def DF(p_c,V_0,T,l_eff,p):         #DERIVEE GREENSHIELFD
        return V_0*(1-2*p*l_eff)
    DF = np.vectorize(DF)
    
    def FINV(V_0,p_max,p):           #inverse de f' de greenshields
        return 0.5*p_max*(1-p/V_0);
    FINV = np.vectorize(FINV)
        
#===================  MODELE GREENSHIELDS 2  ======================================#
if case == 2:
    
    def V(p_c,T,V_0,p_max,l_eff,p,n):  #VITESSE GREENSHIELDS
        return V_0*(1 - (p/p_max)**n_green)
    V = np.vectorize(V)
        
    def F(p_c,T,l_eff,V_0,p):          #FLUX GRRENSHIELDS
        return p*V_0*(1 - (p*l_eff)**n_green)
    F = np.vectorize(F)

    def DF(p_c,V_0,T,l_eff,p):         #DERIVEE GREENSHIELFD
        return V_0*(1 - 3*(p*l_eff)**n_green)
    DF = np.vectorize(DF)

    def FINV(V_0,p_max,p):             #inverse de f' de greenshields
        return p_max*((1-p/V_0)/(n_green +1))**(1/n_green);
    FINV = np.vectorize(FINV)
    
    
#=================== CLASSE Section ================================#
#Un objet de type section correspond à une section de route avec le modele LWR.                       
#CARACTERISTIQUES DE CHAQUES SECTIONS DE ROUTE

class section:
    p_max = None                #Densite maximale veh/km
    V_0 = None                  #Vitesse desiree km/h
    T = None                    #Ecart de temps en heures 
    f = None                    #Fonction flux  vh/h
    v = None                    #Fonction vitesse  km/h
    df = None                   #Fonction derivé du flux  
    I = None                    #Nombre de lignes                                  
    U = None                    #Discretisation de la section  
    flux_entrant = None         #Flux entrant dans la route (calculé grâce à la classe intersection ou donnée par les conditions limites)
    flux_sortant = None         #Flux sortant de la route (calculé grâce à la classe intersection ou donnée par les conditions limites)
    flux_t_entrant=None         #Flux de temps d'attente rentrant dans la route
    sortie = None               #Boolean pour savoir si la route est la derniere (pas de route après)
    feux = None                 #Si la section dispose d'un feux à la fin de la route on met une instance d'un objet de classe feux        
    l_eff = None                #Longeur moyenne entre deux voitures en mode congestion (km)
    p_c = None                  #Densite critique pour laquelle la fonction triangulaire est maximum
    nom = None                  #Nom caracterisant une section
    t_att = None                #Tableau avec les temps d'attente moyen dans chaque maille
    N = None                    #Tableau avec le nombre cumulatif de véhicules passé par chaque maille
    
#Remarque: on peut definir la route de deux façon: ou en fixant la densite critique (p_c) d'où le temps de réaction T sera déduit, ou inversement
#                sect.section(p_max=p_max, V_0=V_0, p_c = p_c2, nom="Albi entree nord 1", I = 2, U = u_init1)
    def __init__(self, p_max, V_0, nom, I, U, flux_entrant=0,flux_t_entrant = 0. , T= None, p_c = None, sortie = False , feux = None):                          
        print("init")      
        print('I =' +str(I))
        print('p_max' + str(p_max))
        self.p_max = I*p_max
        self.V_0 = V_0
        self.l_eff = 1./self.p_max  
        
        if not p_c:
            self.T = T
            self.p_c = 1./(V_0*T+self.l_eff)
        else:
            self.p_c = p_c
            self.T = (1-p_c/self.p_max)/(p_c*V_0)
            print('T = ' + str(self.T))
        if case != 0:
            self.p_c = self.p_max/(n_green+1)**(1/n_green)
            print((n_green+1)**(1/n_green))
        self.nom = nom
        self.I = I
        self.U = np.copy(U)
        self.f = partial(F, self.p_c ,self.T ,self.l_eff,self.V_0)
        self.v = partial(V,self.p_c,self.T,self.V_0,self.p_max,self.l_eff)
        self.df = partial(DF,self.p_c,self.V_0,self.T,self.l_eff)
        if case== 0:
            self.finv = partial(FINV,self.p_c,self.p_max)
        else:
            self.finv = partial(FINV,self.V_0,self.p_max)
        #flux entrant initiale:
        self.flux_entrant = flux_entrant
        #flux sortant initiale nul:
        self.flux_sortant = 0
        self.flux_t_entrant = flux_t_entrant
        self.sortie = sortie
        self.feux = feux
        
        self.t_att = [0 for i in U]
        self.N = [0 for i in U]
        print("La route "+nom+" a p_c = " + str(self.p_c))
        
#=================== CLASSE FEUX ================================#

#Class feux ou vert va etre True si (t % periode) est compris entre t_init et t_final.
#dt est le pas du schema de Godunov et on va incrementer t de dt a chaque calcul (en secondes)
class feux:
    
    def __init__(self, periode, t_init, t_final, dt , t=0):
        self.periode = periode
        self.t_init = t_init
        self.t_final = t_final
        self.dt = dt 
        self.t = t
        if self.t_init==0:
            self.vert = True
        else:
            self.vert = False
    
    def passage(self):
        self.temps = self.t % self.periode
        if self.temps >= self.t_init and self.temps <= self.t_final:
            self.vert = True
            return True
        else:
            self.vert = False
            return False
        
    def changement(self, periode , t_init, t_final):
        self.periode = periode
        self.t_init = t_init
        self.t_final = t_final
    
    
    


#T=1,p_max=150
def V (p):
    return min (50,3600./(1.*p)-3600./(150.))

def V_prime(p):
    if p<3600/74.:
        print('here')
        return 0
    else:
        return -1800./p**2
V_prime = np.vectorize(V_prime)


