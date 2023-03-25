import modele as mod
import sections as s
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import *
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.stats import norm
from scipy.stats import rv_histogram


all_data = [np.random.normal(0, std, size=100) for std in range(1, 4)]

filename = "MC_matrix"

mc_values = np.loadtxt(filename)

fl_fix = 800
init_t = 10
init_x = 1

def density_dist(mc_v, t, x):
    return mc_v[:][x]

axis_color = 'lightgoldenrodyellow'

def distribution():
    m = np.zeros(len(mc_values[:][0]))
    v = np.zeros(len(mc_values[:][0]))
    for i in range(len(mc_values[:][0])):
        m[i] = np.mean(mc_values[:][i])
        v[i] = np.var(mc_values[:][i])
    return m, v

def histogram(mc_v):
    
    M, N = np.shape(mc_v)
    count = np.zeros((N,N))
    bin = np.zeros((N+1,N+1))
    for i in range(len(mc_v[0])):
        c, b = np.histogram(mc_v[:][i], bins=100)
        count[i] = c
        bin[i] = b
        
    return count, bin


def boxplot(mc_v):
    M, N = np.shape(mc_v)
    plt.figure()
    lab = np.linspace(0,N,100)
    plt.boxplot(mc_v, labels = lab)
    #for i in range(N):
       # plt.boxplot(mc_v[:][i])
        
    plt.show()
        
boxplot(mc_values)



counts, bins = histogram(mc_values)

# x_init = 50
# for i in range(0,150,30):
#     plt.figure()
#     print("Bins:" + str(np.shape(bins[i])))
#     plt.hist(bins[i][:-1], bins[i], weights=counts[i])
#     plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)

# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(left=0.25, bottom=0.25)

# Draw the initial plot
# The 'line' variable is used for modifying the line later
line = ax.hist(bins[x_init][:-1], bins[x_init])
ax.set_ylim(0,100)

# Add two sliders for tweaking the parameters

# Draw another slider
x_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axis_color)
x_slider = Slider(x_slider_ax, label='Position [Maille]',valmin=0,valmax=100,valinit=init_x,valstep=1)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(bins[x_init][:-1], bins[x_init], weights = counts[x_init])
    fig.canvas.draw_idle()
x_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    x_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15], facecolor=axis_color)
color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

plt.show()


