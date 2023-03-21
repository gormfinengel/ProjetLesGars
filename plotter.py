import modele as mod
import sections as s
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import *
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.stats import norm
from scipy.stats import rv_histogram

filename = "MC_matrix"

loaded_matrix = np.loadtxt(filename)

M = 100
T = 100
N = 100

original_matrix = np.zeros((M,T,N))

mc_values = loaded_matrix

#mc_values = loaded_matrix.reshape(loaded_matrix.shape[0], loaded_matrix.shape[1] // original_matrix.shape[2], original_matrix.shape[2])

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
    
    counts = np.zeros(np.shape(mc_v))
    
    bins = np.zeros((np.shape(mc_v)[0] + 1, np.shape(mc_v)[0] + 1))
    print(np.shape(bins))
    for i in range(len(mc_v[0])):
        c, b = np.histogram(mc_v[:][i], bins=100)
        counts[i] = c
        bins[i] = b
        
    return counts, bins
        
counts, bins = histogram(mc_values)

x_init = 50

plt.figure()
plt.hist(bins[50][:-1], bins[50], weights=counts[50])
plt.show()

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


