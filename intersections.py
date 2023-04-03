# -*- coding: utf-8 -*-
import numpy as np
from functools import partial
import sections as mod


#=================== FONCTION INTERSECTION ================================#
#Fonction qui prend en arguments les fluxs initiaux (offre et demande), ainsi que la matrice A et le vecteurP , et qui renvoie les fluxs "gamma chapeaux" en intersection
#Parametres: -n: nombre de voie entrantes
#            -m: nombre de voie sortantes
#            -demande: tableau avec flux entrant
#            -offre: tableau avec flux voies sortantes
#            -A : matrice avec pourcentage de passage (A[i,j]=pourcentage de voiture de la route i qui vont a la route j)
#            -P : vecteur avec le droit de passage lorsque l'offre est limitante
#            -verbose: boolean pour afficher des informations dans la console

#Solution:    -resultat: tableaux avec les flux

def flux_intersection(n,m,demande,offre, A, P,temps_att, verbose = 0):  
    resultat = []
    flux_t = [0 for i in range(m)]
    
    #Cas où il y a qu'une voie d'entrée et une voie de sortie
    if (m==1 and n==1):
        resultat.append(min(offre[0],demande[0]))
        resultat.append(min(offre[0],demande[0]))
        flux_t[0] = resultat[1] * temps_att[0]
        return flux_t,resultat
    else:
        D=np.sum(demande)
        O=np.sum(offre)
        
        #Cas où la demande est limitante i.e tout le flux peut passer        
        if O>=D: 
            if verbose != 0:
                print("La demande est limitante \n")
            flux_sortie = np.dot(A,demande)
            A_flux= [[0 for i in A[0]] for j in A]
            recurrence = False
            for i in range(0,m):
                if flux_sortie[i]>offre[i]:
                    recurrence = True
                    route_sature = i
            
        #Cas ou une voie de sortie est saturee et il faut rediriger les voitures vers une autre voie sortante
            if recurrence:
                correction = offre[route_sature]/flux_sortie[route_sature]
                A_nouveau = np.delete(A,route_sature,axis=0)
                
                #Redimensionnement de A_nouveau pour que la somme des colonnes soient égales à 1 puis on rappelle la fonction intersection en enlevant la voie saturée
                redi = [1/(1-j) if (j!=1) else 100 for j in A[route_sature]]
                A_nouveau = [np.array(j)*np.array(redi) for j in A_nouveau]
                m_nouveau = m-1
                offre_nouveau = np.delete(offre,route_sature)
                demande_nouveau = [demande[i]-correction*demande[i]*A[route_sature][i] for i in range(0,n)]
                t,resultat = flux_intersection(n,m_nouveau,demande_nouveau,offre_nouveau,A_nouveau,P, temps_att)
                resultat = np.insert(resultat,n+route_sature,offre[route_sature])
                for i in range(0,n):
                    resultat[i] = resultat[i] + correction*demande[i]*A[route_sature][i]
                if verbose != 0:
                    print(resultat)
                A_flux = np.array(A) * np.array(demande)
                flux_t = [float(np.sum(np.array(i)*np.array(temps_att))) for i in A_flux]
                return flux_t,resultat
            
            #Cas ou aucune voie de sortie est saturee    
            else:
                resultat = np.append(demande,flux_sortie)
                if verbose != 0:
                    print(resultat)
                A_flux = np.array(A) * np.array(demande)
                flux_t = [float(np.sum(np.array(i)*np.array(temps_att))) for i in A_flux]
                return flux_t,resultat
            
        #Cas ou l'offre est trop petite comparé à la demande i.e pas tout le flux passe
        else:
            if verbose != 0:
                print("L'offre est limitante \n")
            flux_entree = np.array(P)*O
            flux_t = [np.dot(A[i],temps_att) for i in range(m)]    
            recurrence = False
            for i in range(0,n):
                if flux_entree[i]>demande[i]:
                    recurrence = True
                    route_vide = i
                
            #Cas ou une voie d'entree à le droit a plus de flux que besoin     
            if (recurrence and n>1):
                A_nouveau = np.delete(A,route_vide,axis=1)
                P_nouveau = np.delete(P,route_vide)
                
                #Redimensionnement de p pour que la somme soit égale à 1 puis on rapelle la fonction intersection en enlevant la route d'entrée déjà asservit.
                if P[route_vide]==1:
                    P_nouveau = [1/(n-1) for i in P_nouveau]
                else:
                    P_nouveau = np.array(P_nouveau)/(1-P[route_vide])
                
                demande_nouveau = np.delete(demande,route_vide)
                offre_nouveau = [offre[i]-demande[route_vide]*A[i][route_vide] for i in range(0,m)]
                n_nouveau=n-1
                t_att_nouveau = np.delete(temps_att,route_vide)
                t,resultat = flux_intersection(n_nouveau,m,demande_nouveau,offre_nouveau,A_nouveau,P_nouveau,t_att_nouveau)
                resultat = np.insert(resultat,route_vide,demande[route_vide])
                for i in range(n,m+n):
                    resultat[i] = resultat[i] + demande[route_vide]*A[i-n][route_vide]
                if verbose != 0:
                    print(resultat)
                flux_t = [flux_t[i] * resultat[n+i] for i in range(m)]
                A_flux = A * np.array(resultat[0:n])
                flux_t = [np.sum(np.array(i)*np.array(temps_att)) for i in A_flux]
                return flux_t,resultat 
          
            #Cas ou aucune voie d'entrée a de droit a plus de flux que besoin             
            else:
                resultat = np.append(flux_entree,offre)
                if verbose != 0:
                    print(resultat)
                flux_t = [flux_t[i] * resultat[n+i] for i in range(m)]
                A_flux = A * np.array(resultat[0:n])
                flux_t = [np.sum(np.array(i)*np.array(temps_att)) for i in A_flux]
                return flux_t,resultat
            
            
#=================== CLASSE INTERSECTION ================================#
#Classe qui modélise une intersection. 

#                ATTRIBUTS
# route_entree:     tableau avec tous les routes en entrée de l'intersection
# route_sortie:     tableau avec tous les routes en sortie de l'intersection
# A:            matrice avec pourcentage de passage (A[i,j]=pourcentage de voiture de la route i qui vont a la route j)
# P:            vecteur avec le droit de passage lorsque l'offre est limitante
# verbose:         boolean pour afficher des informations dans la console
#=========================================================================#
    
            
class intersection:
    route_entree = None
    route_sortie = None
    A = None
    P = None
    verbose = None
    
    def __init__(self,route_entree,route_sortie,A,P,verbose=0):
        self.route_entree = route_entree
        self.route_sortie = route_sortie
        self.A = A
        self.P = P
        self.verbose=verbose
    
    #Calcul le flux d'entree et de sortie grâce à la fonction flux_intersection définit en haut
    def flux(self):
        n = len(self.route_entree)
        m = len(self.route_sortie)
        demande = [mod.demande(i) for i in self.route_entree]
        offre = [mod.offre(i) for i in self.route_sortie]
        t_att = [i.t_att[len(i.U)-1] for i in self.route_entree]
        flux_t,resultat = flux_intersection(n,m,demande,offre,self.A,self.P,t_att, verbose  = self.verbose)
        for i in range(0,n):
            self.route_entree[i].flux_sortant = resultat[i]
        for i in range(0,m):
            self.route_sortie[i].flux_entrant = resultat[i+n]
            self.route_sortie[i].flux_t_entrant = flux_t[i]
        return resultat
    
