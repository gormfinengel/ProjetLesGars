import sections as s
import modele as mod
import intersections as inter
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [15, 5]
plt.rcParams["figure.autolayout"] = True    



def bruit(road, cst_flux, var, scenario=None, iter=None):
    sigma = np.random.normal(0,var)
    if scenario is not None:
        road.flux_entrant = scenario[-(1+iter)] + sigma*500
    else:
        road.flux_entrant = cst_flux+sigma*500

def road_init(N, h, dt, I, p_ref, name, p_max, V_0):
    U_init = np.ones(N)*p_ref
    road = s.section(p_max, V_0, name, I, U_init, sortie=True)
    model = mod.model([road], intersections=None, carrefours=None, dt=dt, h=h)
    return road, model

def simulation(road, model, Niter, iter_tt, scenario, cst_flux, var):

    snapshots = np.linspace(0,Niter-1,5,dtype=int)
    density_stock = []
    numcell_stock = []
    T = []
    iter = 0
    fin_tt = False
    depart_tt = False
    numcell_tt = 1
    vit_tt = 0
    dist_tt = 0
    tocell_tt = 0

    while iter <= Niter or fin_tt == False:
        going_tt = False
        model.resoudre()
        if iter < len(scenario):
            bruit(road, cst_flux, var, scenario, iter)
        else: 
            bruit(road, cst_flux, var)
        if iter == iter_tt:
            depart_tt = True
        if depart_tt and not fin_tt:
            dt = model.dt
            while not going_tt:
                vit_tt = road.V_0*(1-road.U[numcell_tt-1]/road.p_max)
                tocell_tt = (model.h*numcell_tt-dist_tt)/vit_tt
                if dt < tocell_tt:
                    dist_tt += vit_tt*dt
                    going_tt = True
                else: 
                    dist_tt = numcell_tt*model.h
                    numcell_tt += 1
                    dt -= tocell_tt
                if numcell_tt >= len(road.U):
                    fin_tt = True
                    going_tt = True
                    fin_tt_iter = iter
        for s in snapshots:
            if iter==s:
                T.append((iter+1)*model.dt)
                numcell_stock.append(numcell_tt)
                density_stock.append(np.copy(road.U))
        iter += 1
        
    return np.array(density_stock), np.array(T), fin_tt_iter*model.dt, np.array(numcell_stock)

def monte_carlo(M, Niter, N, iter_tt, h, dt, I, p_ref, name, p_max, V_0, scenario, var, cst_flux):
    
    stockage = []
    T = np.zeros(5)
    tt_stockage = []
    numcell_stock = []

    for m in range(M):
        road, model = road_init(N, h, dt, I, p_ref, name, p_max, V_0)
        simu_densities, T, tt_T, numscell = simulation(road, model, Niter, iter_tt, scenario, cst_flux, var)
        tt_stockage.append(tt_T)
        stockage.append(simu_densities)
        numcell_stock.append(numscell)
        print("##################################")
        print("Itération MONTE CARLO: " + str(m))
        print("##################################")

    numscell_mean = np.zeros(5)
    for s in range(5):
        numscell_vals = np.array(numcell_stock)[:,s]
        numscell_mean[s] = np.mean(numscell_vals, axis=0)
        
    return np.asarray(stockage), T, tt_stockage, numscell_mean

def data_org(values):
    
    # Calcul du vecteur moyenne pour chaque snapshot
    m = []
    var = []
    val_to_plot = []
    for s in range(5):
        val = values[:,s]
        m.append(np.mean(val, axis=0))
        var.append(np.std(val, axis=0)**2)
        val_to_plot.append(val)

    return np.array(m), np.array(var), np.array(val_to_plot)

def plot_mc(m, var, values, T, numcell):

    for s in range(np.shape(values)[0]):

        plt.subplot(122)
        plt.title("État du traffic routier à " + str(T[s]*60) + " minutes.")
        sns.boxplot(values[s], showfliers=False)
        plt.scatter(numcell[s], 0, marker = "X", c='r')
        plt.xlabel("Route (en maille)")
        plt.ylabel("Densité (en vh/km)")

        plt.subplot(221)
        plt.xlabel('Route (en maille)', c='r')
        plt.ylabel('Moyenne', c='r')
        plt.plot(m[s], c='r')

        plt.subplot(223)
        plt.xlabel('Route (en maille)', c='b')
        plt.ylabel('Variance', c='b')
        plt.plot(var[s], c='b')


        plt.tight_layout()
        plt.show()

def plot_ref(values, T, numcell):
    for s in range(5):
        plt.title("État du traffic routier à " + str(T[s]*60) + " minutes.")
        plt.plot(values[s,0])
        plt.scatter(numcell[s], 0, marker = "X", c='r')
        plt.show()

def compute_flux(p, V_0, p_max):
    return p*V_0*(1-p/p_max)

def make_sin_flux(lenf, f_r, f_l, nb_per):
    x = np.linspace(0, (1+nb_per*2)*np.pi,lenf)
    if f_r == f_l:
        amplitude = 0.35*f_r
    else:
        amplitude = np.abs(f_r-f_l)/2
    vertical_shift = (f_r+f_l)/2
    phase_shift = np.arcsin((f_r-vertical_shift)/amplitude)
    flux = amplitude*np.sin(x-phase_shift) + vertical_shift
    return flux



#################################################################
######################### SIMULATIONS ###########################
#################################################################

#Route
N = 150
h = .1
V_0 = 80
dt = 0.8*(h/V_0)
I = 1
p_max = 100
name = "Route de test"

var0 = .01
var1 = .15
var2 = .3
var3 = .5
var4 = 1

var = var4

M = 100
N_iter = 300

#FONCTION BRUITEUSE
def bruitage(flux, var, cst_flux):
    sigma = np.random.normal(0, var, len(flux))
    flux_bruite = flux + sigma*500
    return flux_bruite

#Scénario REFERENCE
len_scenar_ref = 100
p_ref_ref = 40
flux_ref_ref = compute_flux(p_ref_ref, V_0, p_max)
cst_flux_ref = flux_ref_ref
iter_tt_ref = 50
scenario_ref = np.ones(len_scenar_ref)*flux_ref_ref

#values_ref, T_ref, tt_stock_ref, numscell_ref = monte_carlo(M, N_iter, N, iter_tt_ref, h, dt, I, p_ref_ref, name, p_max, V_0, scenario_ref, var, cst_flux_ref)
values_ref_ref, T_ref_ref, tt_stock_ref_ref, numscell_ref_ref = monte_carlo(1, N_iter, N, iter_tt_ref, h, dt, I, p_ref_ref, name, p_max, V_0, scenario_ref, 0, cst_flux_ref)
#tt_stock_ref = np.array(tt_stock_ref)*60
#tt_stock_ref_ref = np.array(tt_stock_ref_ref)*60
#tt_mean_ref = np.mean(tt_stock_ref)
#tt_mean_ref_ref = np.mean(tt_stock_ref_ref)
#m_ref, var_ref, values_m_ref = data_org(values_ref)



#Scénario JUMP
len_scenar_jump = 100
p_ref_jump = 12
flux_ref_jump = compute_flux(p_ref_jump, V_0, p_max)
cst_flux_jump = 2.*flux_ref_jump
iter_tt_jump = 50
scenario_jump = np.multiply(np.ones(len_scenar_jump)*flux_ref_jump, np.linspace(2.,1.,len_scenar_jump))

#values_jump, T_jump, tt_stock_jump, numscell_jump = monte_carlo(M, N_iter, N, iter_tt_jump, h, dt, I, p_ref_jump, name, p_max, V_0, scenario_jump, var, cst_flux_jump)
values_jump_ref, T_jump_ref, tt_stock_jump_ref, numscell_jump_ref = monte_carlo(1, N_iter, N, iter_tt_jump, h, dt, I, p_ref_jump, name, p_max, V_0, scenario_jump, 0, cst_flux_jump)
#tt_stock_jump = np.array(tt_stock_jump)*60
#tt_stock_jump_ref = np.array(tt_stock_jump_ref)*60
#tt_mean_jump = np.mean(tt_stock_jump)
#tt_mean_jump_ref = np.mean(tt_stock_jump_ref)
#m_jump, var_jump, values_m_jump = data_org(values_jump)


#Scénario DROP
len_scenar_drop = 100
p_ref_drop = 30
flux_ref_drop = compute_flux(p_ref_drop, V_0, p_max)
cst_flux_drop = .4*flux_ref_drop
iter_tt_drop = 50
scenario_drop = np.multiply(np.ones(len_scenar_drop)*flux_ref_drop, np.linspace(.4,1.,len_scenar_drop))

#values_drop, T_drop, tt_stock_drop, numscell_drop = monte_carlo(M, N_iter, N, iter_tt_drop, h, dt, I, p_ref_drop, name, p_max, V_0, scenario_drop, var, cst_flux_drop)
values_drop_ref, T_drop_ref, tt_stock_drop_ref, numscell_drop_ref = monte_carlo(1, N_iter, N, iter_tt_drop, h, dt, I, p_ref_drop, name, p_max, V_0, scenario_drop, 0, cst_flux_drop)
#tt_stock_drop = np.array(tt_stock_drop)*60
#tt_stock_drop_ref = np.array(tt_stock_drop_ref)*60
#tt_mean_drop = np.mean(tt_stock_drop)
#tt_mean_drop_ref = np.mean(tt_stock_drop_ref)
#m_drop, var_drop, values_m_drop = data_org(values_drop)



#Scénario SINUS JUMP
len_scenar_sinus_jump = 100
nb_per_jump = 1
p_ref_sin_jump = 12
f_r_jump = compute_flux(p_ref_sin_jump, V_0, p_max)
f_l_jump = 2*f_r_jump
iter_tt_sin_jump = 50
scenario_sin_jump = make_sin_flux(len_scenar_sinus_jump, f_r_jump, f_l_jump, nb_per_jump)

#values_sin_jump, T_sin_jump, tt_stock_sin_jump, numscell_sin_jump = monte_carlo(M, N_iter, N, iter_tt_sin_jump, h, dt, I, p_ref_sin_jump, name, p_max, V_0, scenario_sin_jump, var, f_l_jump)
values_sin_jump_ref, T_sin_jump_ref, tt_stock_sin_jump_ref, numscell_sin_jump_ref = monte_carlo(1, N_iter, N, iter_tt_sin_jump, h, dt, I, p_ref_sin_jump, name, p_max, V_0, scenario_sin_jump, 0, f_l_jump)
#tt_stock_sin_jump = np.array(tt_stock_sin_jump)*60
#tt_stock_sin_jump_ref = np.array(tt_stock_sin_jump_ref)*60
#tt_mean_sin_jump = np.mean(tt_stock_sin_jump)
#tt_mean_sin_jump_ref = np.mean(tt_stock_sin_jump_ref)
#m_sin_jump, var_sin_jump, values_m_sin_jump = data_org(values_sin_jump)



#Scénario SINUS DROP
len_scenar_sinus_drop = 100
nb_per_drop = 1
p_ref_sin_drop = 30
f_r_drop = compute_flux(p_ref_sin_drop, V_0, p_max)
f_l_drop = .4*f_r_drop
iter_tt_sin_drop = 50
scenario_sin_drop = make_sin_flux(len_scenar_sinus_drop, f_r_drop, f_l_drop, nb_per_drop)

#values_sin_drop, T_sin_drop, tt_stock_sin_drop, numscell_sin_drop = monte_carlo(M, N_iter, N, iter_tt_sin_drop, h, dt, I, p_ref_sin_drop, name, p_max, V_0, scenario_sin_drop, var, f_l_drop)
values_sin_drop_ref, T_sin_drop_ref, tt_stock_sin_drop_ref, numscell_sin_drop_ref = monte_carlo(1, N_iter, N, iter_tt_sin_drop, h, dt, I, p_ref_sin_drop, name, p_max, V_0, scenario_sin_drop, 0, f_l_drop)
#tt_stock_sin_drop = np.array(tt_stock_sin_drop)*60
#tt_stock_sin_drop_ref = np.array(tt_stock_sin_drop_ref)*60
#tt_mean_sin_drop = np.mean(tt_stock_sin_drop)
#tt_mean_sin_drop_ref = np.mean(tt_stock_sin_drop_ref)
#m_sin_drop, var_sin_drop, values_m_sin_drop = data_org(values_sin_drop)


#Scénario SINUS REFERENCE
len_scenar_sinus_ref = 100
nb_per_ref = 1
p_ref_sin_ref = 30
f_r_ref = compute_flux(p_ref_sin_ref, V_0, p_max)
f_l_ref = f_r_ref
iter_tt_sin_ref = 50
scenario_sin_ref = make_sin_flux(len_scenar_sinus_ref, f_r_ref, f_l_ref, nb_per_ref)

#values_sin_ref, T_sin_ref, tt_stock_sin_ref, numscell_sin_ref = monte_carlo(M, N_iter, N, iter_tt_sin_ref, h, dt, I, p_ref_sin_ref, name, p_max, V_0, scenario_sin_ref, var, f_l_ref)
values_sin_ref_ref, T_sin_ref_ref, tt_stock_sin_ref_ref, numscell_sin_ref_ref = monte_carlo(1, N_iter, N, iter_tt_sin_ref, h, dt, I, p_ref_sin_ref, name, p_max, V_0, scenario_sin_ref, 0, f_l_ref)
#tt_stock_sin_ref = np.array(tt_stock_sin_ref)*60
#tt_stock_sin_ref_ref = np.array(tt_stock_sin_ref_ref)*60
#tt_mean_sin_ref = np.mean(tt_stock_sin_ref)
#tt_mean_sin_ref_ref = np.mean(tt_stock_sin_ref_ref)
#m_sin_ref, var_sin_ref, values_m_sin_ref = data_org(values_sin_ref)

fig, axs = plt.subplots(2,3, figsize=(15,7))

plt.suptitle("Density outputs for non-noisy flux at T = 13.5 minutes")

axs[0,0].plot(values_ref_ref[0,3])
axs[0,0].set_ylim(30,50)
axs[0,0].set_xlabel("Road (in meshes)")
axs[0,0].set_ylabel("Density (in vh/km)")
axs[0,0].set_title("REFERENCE output")

axs[0,1].plot(values_jump_ref[0,3])
axs[0,1].set_xlabel("Road (in meshes)")
axs[0,1].set_ylabel("Density (in vh/km)")
axs[0,1].set_title("JUMP output")

axs[0,2].plot(values_drop_ref[0,3])
axs[0,2].set_xlabel("Road (in meshes)")
axs[0,2].set_ylabel("Density (in vh/km)")
axs[0,2].set_title("DROP output")

axs[1,0].plot(values_sin_ref_ref[0,3])
axs[1,0].set_xlabel("Road (in meshes)")
axs[1,0].set_ylabel("Density (in vh/km)")
axs[1,0].set_title("SIN REFERENCE output")

axs[1,1].plot(values_sin_jump_ref[0,3])
axs[1,1].set_xlabel("Road (in meshes)")
axs[1,1].set_ylabel("Density (in vh/km)")
axs[1,1].set_title("SIN JUMP output")

axs[1,2].plot(values_sin_drop_ref[0,3])
axs[1,2].set_xlabel("Road (in meshes)")
axs[1,2].set_ylabel("Density (in vh/km)")
axs[1,2].set_title("SIN DROP output")

plt.show()


###############################################################
########################## TT #################################
###############################################################

'''
TT = np.array([tt_stock_ref, tt_stock_jump, tt_stock_drop, tt_stock_sin_ref, tt_stock_sin_jump, tt_stock_sin_drop])

plt.boxplot(TT.T)
plt.ylabel("travel time (in minutes)")
plt.xticks([1,2,3,4,5,6],["REFERENCE", "JUMP", "DROP", "SIN REFERENCE", "SIN JUMP", "SIN DROP"])
plt.xlabel("Flux scenarios")
plt.title("Travel time distributions with variance = " + str(var))
plt.show()

print("Temps de trajet de référence pour un véhicule (flux non-bruité): " + str(tt_mean_ref_ref))
print("Temps de trajet de référence pour un véhicule (flux bruité): " + str(tt_mean_ref))
print("Variance temps de trajet de référence pour un véhicule (flux bruité): " + str(np.std(tt_stock_ref)**2))

print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux non-bruité): " + str(tt_mean_jump_ref))
print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux bruité): " + str(tt_mean_jump))
print("Variance temps de trajet pour un véhicule avec un départ après le JUMP (flux bruité): " + str(np.std(tt_stock_jump)**2))

print("Temps de trajet pour un véhicule avec un départ après le DROP (flux non-bruité): " + str(tt_mean_drop_ref))
print("Temps de trajet pour un véhicule avec un départ après le DROP (flux bruité): " + str(tt_mean_drop))
print("Variance temps de trajet pour un véhicule avec un départ après le DROP (flux bruité): " + str(np.std(tt_stock_drop)**2))

print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux sinus non-bruité): " + str(tt_mean_sin_jump_ref))
print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux sinus bruité): " + str(tt_mean_sin_jump))
print("Variance temps de trajet pour un véhicule avec un départ après le JUMP (flux sinus bruité): " + str(np.std(tt_stock_sin_jump)**2))

print("Temps de trajet pour un véhicule avec un départ après le DROP (flux sinus non-bruité): " + str(tt_mean_sin_drop_ref))
print("Temps de trajet pour un véhicule avec un départ après le DROP (flux sinus bruité): " + str(tt_mean_sin_drop))
print("Variance temps de trajet pour un véhicule avec un départ après le DROP (flux sinus bruité): " + str(np.std(tt_stock_sin_drop)**2))

print("Temps de trajet de référence pour un véhicule (flux sinus non-bruité): " + str(tt_mean_sin_ref_ref))
print("Temps de trajet de référence pour un véhicule (flux sinus bruité): " + str(tt_mean_sin_ref))
print("Variance temps de trajet de référence pour un véhicule (flux sinus bruité): " + str(np.std(tt_stock_sin_ref)**2))



###############################################################
######################### PLOTS ###############################
###############################################################

fig, axs = plt.subplots(2,3,figsize=(15,8))

axs[0,0].plot(bruitage(scenario_ref, var, flux_ref_ref), c='b')
axs[0,0].set_ylim(1300, 2800)
axs[0,0].set_title("Reference flux with noise")

axs[0,1].plot(bruitage(scenario_jump, var, flux_ref_jump), c='r')
axs[0,1].set_title("Jump flux with noise")

axs[0,2].plot(bruitage(scenario_drop, var, cst_flux_drop), c='orange')
axs[0,2].set_title("Drop flux with noise")

axs[1,0].plot(bruitage(scenario_sin_ref, var, f_r_ref), c='green')
axs[1,0].set_title("Sinusoidal reference flux with noise")

axs[1,1].plot(bruitage(scenario_sin_drop, var, f_l_drop), c='black')
axs[1,1].set_title("Sinusoidal drop flux with noise")

axs[1,2].plot(bruitage(scenario_sin_jump, var, f_r_jump), c='pink')
axs[1,2].set_title("Sinusoidal jump flux with noise")

for i in range(2):
    for j in range(3):
        axs[i,j].set_xlabel("Road (in cells of size h)")
        axs[i,j].set_ylabel("Flux (in vehicles per hour)")

plt.tight_layout()
plt.show()


###################################################################

for s in range(5):

    fig, axs = plt.subplots(2,3,figsize=(15,10))

    fig.suptitle("Uncertainties propagation at " + str(T_ref[s]*60) + " minutes, with variance = " +str(var))

    axs[0,0].boxplot(values_m_ref[s])
    axs[0,0].plot(values_ref_ref[0][s], label="Reference density", c='g')
    axs[0,0].scatter(numscell_ref[s], 0, marker = "X", c='r')
    axs[0,0].scatter(numscell_ref_ref[s], 0, marker = "X", c='g')
    axs[0,0].set_title("Reference flux output")
    axs[0,0].set_xticks([0,20,40,60,80,100,120,140],[0,20,40,60,80,100,120,140])
    axs[0,0].legend(loc=0)
    #axs00twin = axs[0,0].twinx()
    #axs00twin.plot(var_ref[s], label="density variance")
    #axs00twin.legend(loc=7)


    axs[0,1].boxplot(values_m_jump[s])
    axs[0,1].plot(values_jump_ref[0][s], label="Reference density", c='g')
    axs[0,1].scatter(numscell_jump[s], 0, marker = "X", c='r')
    axs[0,1].scatter(numscell_jump_ref[s], 0, marker = "X", c='g')
    axs[0,1].legend(loc=0)
    axs[0,1].set_xticks([0,20,40,60,80,100,120,140],[0,20,40,60,80,100,120,140])
    axs[0,1].set_title("JUMP flux output")
    #axs01twin = axs[0,1].twinx()
    #axs01twin.plot(var_jump[s], label="density variance")
    #axs01twin.legend(loc=7)    


    axs[0,2].boxplot(values_m_drop[s])
    axs[0,2].plot(values_drop_ref[0][s], label="Reference density", c='g')
    axs[0,2].scatter(numscell_drop[s], 0, marker = "X", c='r')
    axs[0,2].scatter(numscell_drop_ref[s], 0, marker = "X", c='g')
    axs[0,2].set_title("DROP flux output")
    axs[0,2].set_xticks([0,20,40,60,80,100,120,140],[0,20,40,60,80,100,120,140])
    axs[0,2].legend(loc=0)
    #axs02twin = axs[0,2].twinx()
    #axs02twin.plot(var_drop[s], label="density variance")
    #axs02twin.legend(loc=7)
    

    axs[1,0].boxplot(values_m_sin_ref[s])
    axs[1,0].plot(values_sin_ref_ref[0][s], label="Reference density", c='g')
    axs[1,0].scatter(numscell_sin_ref[s], 0, marker = "X", c='r')
    axs[1,0].scatter(numscell_sin_ref_ref[s], 0, marker = "X", c='g')
    axs[1,0].set_title("Sinusoidal reference flux output")
    axs[1,0].set_xticks([0,20,40,60,80,100,120,140],[0,20,40,60,80,100,120,140])
    axs[1,0].legend(loc=0)
    #axs10twin = axs[1,0].twinx()
    #axs10twin.plot(var_sin_ref[s], label="density variance")
    #axs10twin.legend(loc=7)
    

    axs[1,1].boxplot(values_m_sin_jump[s])
    axs[1,1].plot(values_sin_jump_ref[0][s], label="Reference density", c='g')
    axs[1,1].scatter(numscell_sin_jump[s], 0, marker = "X", c='r')
    axs[1,1].scatter(numscell_sin_jump_ref[s], 0, marker = "X", c='g')
    axs[1,1].set_title("Sinusoidal JUMP flux output")
    axs[1,1].set_xticks([0,20,40,60,80,100,120,140],[0,20,40,60,80,100,120,140])
    axs[1,1].legend(loc=0)
    #axs11twin = axs[1,1].twinx()
    #axs11twin.plot(var_sin_jump[s], label="density variance")
    #axs11twin.legend(loc=7)
    

    axs[1,2].boxplot(values_m_sin_drop[s])
    axs[1,2].plot(values_sin_drop_ref[0][s], label="Reference density", c='g')
    axs[1,2].scatter(numscell_sin_drop[s], 0, marker = "X", c='r')
    axs[1,2].scatter(numscell_sin_drop_ref[s], 0, marker = "X", c='g')
    axs[1,2].set_title("Sinusoidal DROP flux output")
    axs[1,2].set_xticks([0,20,40,60,80,100,120,140],[0,20,40,60,80,100,120,140])
    axs[1,2].legend(loc=0)
    #axs12twin = axs[1,2].twinx()
    #axs12twin.plot(var_sin_drop[s], label="density variance")
    #axs12twin.legend(loc=7)
    


    plt.tight_layout()
    plt.show()

###################################################################
######################## TRAVEL TMIE ##############################
###################################################################

print("Temps de trajet de référence pour un véhicule (flux non-bruité): " + str(tt_mean_ref_ref))
print("Temps de trajet de référence pour un véhicule (flux bruité): " + str(tt_mean_ref))

print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux non-bruité): " + str(tt_mean_jump_ref))
print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux bruité): " + str(tt_mean_jump))

print("Temps de trajet pour un véhicule avec un départ après le DROP (flux non-bruité): " + str(tt_mean_drop_ref))
print("Temps de trajet pour un véhicule avec un départ après le DROP (flux bruité): " + str(tt_mean_drop))

print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux sinus non-bruité): " + str(tt_mean_sin_jump_ref))
print("Temps de trajet pour un véhicule avec un départ après le JUMP (flux sinus bruité): " + str(tt_mean_sin_jump))

print("Temps de trajet pour un véhicule avec un départ après le DROP (flux sinus non-bruité): " + str(tt_mean_sin_drop_ref))
print("Temps de trajet pour un véhicule avec un départ après le DROP (flux sinus bruité): " + str(tt_mean_sin_drop))

print("Temps de trajet de référence pour un véhicule (flux sinus non-bruité): " + str(tt_mean_sin_ref_ref))
print("Temps de trajet de référence pour un véhicule (flux sinus bruité): " + str(tt_mean_sin_ref))
'''
