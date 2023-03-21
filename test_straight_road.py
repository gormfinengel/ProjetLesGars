import modele as mod
import sections as s
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import *
from matplotlib.widgets import Slider, Button, RadioButtons

M = 100
T = 100
N = 100

fl_fix = 800

#Variation du flux d'entrée
def variation(route, fl_fix):
    sigma = np.random.normal(0,1)*100
    route.flux_entrant = fl_fix + sigma
    # sigma = np.random.normal(0,0.1)
    # if (int(route.flux_entrant + sigma*fl_fix)>0 and int(route.flux_entrant + sigma*fl_fix)<1600):
    #     route.flux_entrant = route.flux_entrant + sigma*fl_fix
    
               
def init_route():
    U_init = np.zeros(N)           #Discrétisation de la route
    #U_init[0:30] = 100
    h = 0.1                         #Taille de chaque maille (en km)
    I = 1                          #Nombre de lignes
    nom = "Section solo test"
    p_c = 100.                     #Densité critique (en vh/km)
    p_max = 150.                   #Densité max (en vh/km)
    V_0 = 50.                      #Vitesse souhaitée (en km/h)                         
    flux_entrant = 800            #Flux de véhicules entrant (en vh/h)
    fl_fix = flux_entrant
    #flux_t_entrant = 0            #Flux de temps d'attente
    #T = 1.4/3600                   #Temps de réaction
    dt =  0.99*(h/V_0)             #CFL

    route = s.section(p_max = p_max, V_0 = V_0, I = I, U = U_init, flux_entrant=flux_entrant, p_c = p_c, nom = nom, sortie=True)
    model = mod.model([route], intersections= None, carrefours= None, dt = dt, h = h)

    return route, model

def simulation(T):
    # T correspond au nombre de pas de temps qu'on veut pour notre simulation
    
    route, model = init_route()
    
    for i in range(T):
        model.resoudre()
        variation(route,fl_fix)
        
    return route.U

# def simulation(T):
#     # T correspond au nombre de pas de temps qu'on veut pour notre simulation
#     pho_global = np.zeros((T-100,N))
#     route, model = init_route()

#     for i in range(T):
#         model.resoudre()
#         variation(route, fl_fix)
#         if i>=100:
#             pho_global[i-100] = route.U

#     return np.asarray(pho_global)

def MC_sim(M, simulator, T, N):
    # M correspond au nombre d'iter pour MC
    # simulator correspond au modèle à simuler
    stockage = []
    for m in range(M):
        sim = simulator(T)
        stockage.append(sim)
    print(np.shape(stockage))
    return np.asarray(stockage)

mc_values = MC_sim(M, simulation, T, N)

#mc_values_reshaped = mc_values.reshape(mc_values.shape[0], -1)

np.savetxt("MC_matrix",mc_values)

filename = "MC_matrix"

