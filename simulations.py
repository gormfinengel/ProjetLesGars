import sections as s
import modele as mod
import intersections as inter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def following_roads(Nb_R, N, h, I, names, p_c, p_max, V_0, input_fl, exits):

    roads = []

    for r in range(Nb_R):
        U_init = np.zeros(N[r])
        if not exits[r]:
            roads.append(s.section(p_max = p_max[r], V_0 = V_0[r], I = I[r], U = U_init, flux_entrant=input_fl[0], p_c = p_c[r], nom = names[r]))
        elif exits[r]:
            roads.append(s.section(p_max = p_max, V_0 = V_0[r], I = I[r], U = U_init, nom = names[r], p_c = 50., sortie = True))


    #intersection = inter.intersection([roads[0]], [roads[1]], A=[[1.]], P=[1])


    dt = .99*(h/V_0[0]) 
    model = mod.model(roads, intersections=None, carrefours= None, dt = dt, h = h)

    return roads, model

def single_road(N, h, I, names, p_c, p_max, V_0, input_fl, exits):
    
    roads = []

    U_init = np.zeros(N[0])
    roads.append(s.section(p_max = p_max[0], V_0 = V_0[0], I = I[0], U = U_init, flux_entrant=input_fl[0], p_c = p_c[0], nom = names[0], sortie = exits[0]))

    dt = .99*(h/V_0[0]) 
    model = mod.model(roads, intersections=None, carrefours= None, dt = dt, h = h)

    return roads, model


def bruit(road, ind, fl_vec, cst_fl,var):
    sigma = np.random.normal(0,var)
    road.flux_entrant = fl_vec[ind] + sigma*cst_fl

def fl_shock(beg_fl, end_fl, Niter, *, pow = 1):

    input_flow = np.zeros(Niter)

    if Niter/(pow*10) <= Niter:
        trans = Niter/(pow*10)
    else :
        trans = Niter

    dir_c = (end_fl-beg_fl)/trans

    x = 0
    for i in range(Niter):
        if i < (Niter-trans)/2:
            input_flow[i] = beg_fl
        elif i >= (Niter-trans)/2 and i <= (Niter+trans)/2:
            input_flow[i] = x*dir_c + beg_fl
            x+=1
        elif i > (Niter+trans)/2:
            input_flow[i] = end_fl
    
    return input_flow


def fl_cb(beg_fl, end_fl, mid_fl, Niter, rp, *, plateau = False, plateau_size=1):

    input_flow = np.zeros(Niter)

    x1 = 0
    x2 = 0

    if not plateau:
        
        dir_1 = (mid_fl-beg_fl)/rp
        dir_2 = (end_fl-mid_fl)/(Niter-rp)

        for i in range(Niter):
            if i < rp:
                input_flow[i] = x1*dir_1 + beg_fl
                x1 += 1
            elif i == rp:
                input_flow[i] = mid_fl
            elif i > rp:
                input_flow[i] = x2*dir_2 + mid_fl
                x2 += 1
    
    elif plateau:

        if plateau_size==1:
            plateau_size  = Niter/20

        dir_1 = (mid_fl-beg_fl)/(rp-plateau_size/2)
        dir_2 = (end_fl-mid_fl)/(Niter-(rp+plateau_size/2))

        for i in range(Niter):
            if i < rp-plateau_size/2:
                input_flow[i] = x1*dir_1 + beg_fl
                x1 += 1
            elif i >= rp-plateau_size/2 and i <= rp+plateau_size/2:
                input_flow[i] = mid_fl
            elif i > rp+plateau_size/2:
                input_flow[i] = x2*dir_2 + mid_fl
                x2 += 1
    
    return input_flow

def simulation(roads, model, input_fl, cst_fl, var):

    Niter = len(input_fl)

    for i in range(Niter):
        model.resoudre()
        bruit(roads[0], i, input_fl, cst_fl, var)

    densities = []
    for r in roads:
        densities.append(r.U)
    
    return densities


def monte_carlo(M, roads, model, input_fl, cst_fl, var):

    stockage = []
    length = 0
    lengths = []
    for r in roads:
            length += len(r.U)
            lengths.append(len(r.U))

    for m in range(M):
        simu_densities = simulation(roads, model, input_fl, cst_fl, var)
        density_vec = np.zeros(length)
        old_l = 0
        n = 0
        for l in lengths:
            density_vec[old_l:old_l+l] = simu_densities[n]
            old_l = l
            n += 1
        stockage.append(density_vec)

    return np.asarray(stockage)


def plot_mc(values):

    mu = np.mean(values, axis=0)
    var = np.std(values, axis=0)**2

    fig = plt.figure(figsize=(15,10))

    plt.subplot(122)
    box_val = np.zeros((np.shape(values)[1],np.shape(values)[0]))
    for i in range(len(values[0])):
        box_val[i,:] = values[:,i]
    plt.boxplot(box_val.tolist())
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

input_fl_sin = np.linspace(0,25,200)
input_fl_sin = np.sin(input_fl_sin)*200+1700

plt.plot(input_fl_sin)
plt.show()

Nb_R = 1
N = [100]
h = .1
I = [1]
names = ["just a road"]
p_c = [100]
p_max = [150]
V_0 = [50]
exits = [True]

#ref_roads, ref_model = single_road(N=N, h=h, I=I, names=names, p_c=p_c, p_max=p_max, V_0=V_0, input_fl=input_fl_ref, exits=exits)
#plot_mc(monte_carlo(1000, ref_roads, ref_model, input_fl_ref, 800, .05))

#ths_roads, ths_model = single_road(N=N, h=h, I=I, names=names, p_c=p_c, p_max=p_max, V_0=V_0, input_fl=input_fl_ths, exits=exits)
#plot_mc(monte_carlo(1000, ths_roads, ths_model, input_fl_ref, 800, .1))

#sin_roads, sin_model = single_road(N=N, h=h, I=I, names=names, p_c=p_c, p_max=p_max, V_0=V_0, input_fl=input_fl_sin, exits=exits)
#plot_mc(monte_carlo(100, sin_roads, sin_model, input_fl_sin, 1000, .1))


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



        
        """