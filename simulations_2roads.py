from pickle import FALSE
import sections as s
import modele as mod
import intersections as inter
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True    

def noise(flux, p_cst, var):
    flux_noised = np.copy(flux)
    for i in range(len(flux)):
        sigma = np.random.normal(0,var)
        flux_noised[i] = flux[i] + p_cst*sigma
    return flux_noised


def two_roads(N1, N2, h, I, p_ref_1, p_ref_2, name1, name2, p_max_1, p_max_2, V_0, flux_scenario, var):

    roads = []

    U_init_1 = np.ones(N1)*p_ref_1

    U_init_2 = np.ones(N2)*p_ref_2
    flux_scenario_noised = noise(flux_scenario, p_ref_1, var)
    U_init_2[N2-len(flux_scenario):] = flux_scenario_noised

    roads.append(s.section(p_max=p_max_2, V_0=V_0, I=I, U=U_init_2, nom=name2, sortie=False))
    roads.append(s.section(p_max=p_max_1, V_0=V_0, I=I, U=U_init_1, nom=name1, sortie=True))
    intersection = inter.intersection([roads[0]], [roads[1]], A=[[1.]], P=[1])
    
    dt = .99*(h/V_0) 
    model = mod.model(roads, intersections=[intersection], carrefours= None, dt = dt, h = h)

    return roads, model

def simulation(roads, model, Niter):
     
    density_end_simu = roads[0].U[0] 
    print(density_end_simu)
    itermax = Niter
    #if Niter > len(roads[1].U):
       #itermax = len(roads[1].U)

    snapshots = np.linspace(0,itermax-1,5,dtype=int)
    density_stock = []

    for i in range(itermax):
        model.resoudre()
        roads[0].flux_entrant = roads[0].V_0*(1-density_end_simu/roads[0].p_max)
        for s in snapshots:
            if i==s:
                dens_vec = np.zeros(len(roads[0].U)+len(roads[1].U))
                dens_vec[:len(roads[0].U)] = roads[0].U
                dens_vec[len(roads[0].U):] = roads[1].U
                density_stock.append(dens_vec)

    return np.array(density_stock)


def monte_carlo(M, Niter, N1, N2, h, I, p_ref_1, p_ref_2, name1, name2, p_max_1, p_max_2, V_0, flux_scenario, var):
    
    stockage = []

    for m in range(M):
        roads, model = two_roads(N1, N2, h, I, p_ref_1, p_ref_2, name1, name2, p_max_1, p_max_2, V_0, flux_scenario, var)
        simu_densities = simulation(roads, model, Niter)
        
        stockage.append(simu_densities)
        
    return np.asarray(stockage)

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

    return (np.array(m), np.array(var), np.array(val_to_plot))
        

def plot_mc(m, var, values, sep):

    for s in range(np.shape(values)[0]):

        plt.subplot(122)
        sns.boxplot(values[s])
        plt.axvline(x = sep, color = 'b', label = 'Intersection')
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

def plot_ref(values):
    for s in range(5):
        plt.plot(values[s,0])
        plt.show()



scenario_test = np.ones(50)*100
m, var, values = data_org(monte_carlo(10, 300, 100, 100, 0.1, 1, 30, 30, "Route étudiée", "Route expérience", 75., 75., 50., scenario_test, 0.5))
plot_mc(m,var,values,50)

'''


def plot_mc(values):
    
    mu = np.mean(values, axis=0)
    var = np.std(values, axis=0)**2

    fig = plt.figure(figsize=(15,10))

    plt.subplot(122)
    box_val = np.zeros((np.shape(values)[1],np.shape(values)[0]))
    for i in range(len(values[0])):
        box_val[i,:] = values[:,i]
    #plt.plot(box_val.tolist())
    sns.boxplot(box_val.tolist(), showfliers=False)
    plt.xticks([0,20,40,60,80,100])
    plt.xlabel("Route (en maille)")
    plt.ylabel("Densité (en vh/km)")

    plt.subplot(221)
    color = 'tab:red'
    plt.xlabel('Route (en maille)')
    plt.ylabel('Moyenne', color=color)
    plt.plot(mu, color=color)

    plt.subplot(223)
    color = 'tab:blue'
    plt.xlabel('Route (en maille)')
    plt.ylabel('Variance', color=color)
    plt.plot(var, color=color)


    fig.tight_layout()
    plt.show()

####################################################
####################################################


########### REFERENCE SIMULATION ##############

#input_fl_ref = np.ones(200)
#input_fl_ref = input_fl_ref*800

#input_fl_ths = np.ones(200)
#input_fl_ths = input_fl_ref*1850

frequence = 50
temps = np.linspace(0,50)

input_fl_sin = np.sin(2*np.pi*temps*frequence)*500+1300

plt.plot(input_fl_sin)
plt.xlabel("Distance (en maille)")
plt.ylabel("Flux (en vh/min)")
plt.title("Flux entrant ")
plt.show()

Nb_R = 1
N = [100]
h = .1
I = [1]
names = ["just a road"]
p_c = [50]
p_max = [100]
V_0 = [50]
exits = [True]

#ref_roads, ref_model = single_road(N=N, h=h, I=I, names=names, p_c=p_c, p_max=p_max, V_0=V_0, input_fl=input_fl_ref, exits=exits)
#plot_mc(monte_carlo(1000, ref_roads, ref_model, input_fl_ref, 800, .05))

#ths_roads, ths_model = single_road(N=N, h=h, I=I, names=names, p_c=p_c, p_max=p_max, V_0=V_0, input_fl=input_fl_ths, exits=exits)
#plot_mc(monte_carlo(1000, ths_roads, ths_model, input_fl_ref, 800, .1))

values, flux = monte_carlo(10, input_fl_sin, 500, 0, N=N, h=h, I=I, p_ref=50, names=names, p_c=p_c, p_max=p_max, V_0=V_0, exits=exits)
plot_mc(values, flux)


""" 

fig = plt.figure(figsize=(10,8))
ax = plt.axes(xlim=(0, 150), ylim=(0, 150))
line, = plt.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):

    test_model.resoudre()
    variation(test_roads[0], i, input_fl, 800, .5)

    densities = []

    for r in test_roads:
        densities.append(r.U)
    
    plot_densities = np.zeros(150)
    plot_densities[:50] = densities[0]
    plot_densities[50:150] = densities[1]

    x = np.arange(0,150,1)
    y = plot_densities

    line.set_data(x, y)
    
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = FuncAnimation(fig, animate, init_func=init, frames=200, interval=10, blit=True)

anim.save('densite_evolution.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()



'''