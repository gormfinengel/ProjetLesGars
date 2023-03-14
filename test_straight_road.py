import modele as mod
import sections as s
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import *
from matplotlib.widgets import Slider, Button, RadioButtons

M = 100
T = 200
N = 100

fl_fix = 800

init_t = 10
init_x = 49

#Variation du flux d'entrée
def variation(route, fl_fix):
    sigma = np.random.normal(0,0.1)
    if (int(route.flux_entrant + sigma*fl_fix)>0 and int(route.flux_entrant + sigma*fl_fix)<1600):
        route.flux_entrant = route.flux_entrant + sigma*fl_fix
    
               
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
    pho_global = np.zeros((T-N,N))
    route, model = init_route()

    for i in range(T):
        model.resoudre()
        variation(route, fl_fix)
        if i>=100:
            pho_global[i-100] = route.U

    return np.asarray(pho_global)

def MC_sim(M, simulator, T, N):
    # M correspond au nombre d'iter pour MC
    # simulator correspond au modèle à simuler
    stockage = []
    for m in range(M):
        sim = simulator(T)
        stockage.append(sim)
    return stockage

mc_values = MC_sim(M, simulation, T, N)


def density_dist(mc_v, t, x):
    return mc_v[:][t][x]

axis_color = 'lightgoldenrodyellow'

fig = plt.figure()
ax = fig.add_subplot(111)

# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(left=0.25, bottom=0.25)


# Draw the initial plot
# The 'line' variable is used for modifying the line later
[line] = ax.plot(density_dist(mc_values, init_t, init_x), linewidth=2, color='red')
ax.set_ylim(0,100)

# Add two sliders for tweaking the parameters

# Define an axes area and draw a slider in it
t_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
t_slider = Slider(t_slider_ax, label='Time [dt]',valmin=0,valmax=100,valinit=init_t, valstep = 1)

# Draw another slider
x_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axis_color)
x_slider = Slider(x_slider_ax, label='Position [Maille]',valmin=0,valmax=100,valinit=init_x,valstep=1)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(density_dist(mc_values,t_slider.val, x_slider.val))
    fig.canvas.draw_idle()
t_slider.on_changed(sliders_on_changed)
x_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    x_slider.reset()
    t_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15], facecolor=axis_color)
color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

plt.show()


