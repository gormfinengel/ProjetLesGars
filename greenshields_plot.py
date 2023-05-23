import matplotlib.pyplot as plt
import numpy as np

def flow(p, p_max, V_0, n_green):
    return V_0*p*(1-(p/p_max))

def plot_greenshields(p_max, V_0, n_green):

    fig, axs = plt.subplots(1,2,figsize=(15,5))

    p = np.arange(0, p_max)
    p_c = (p_max)/(n_green+1)**(1/n_green)
    f = flow(p, p_max, V_0, n_green)
    f_max = flow(p_c, p_max, V_0, n_green)
    axs[0].plot(p,f)
    axs[0].set_ylim(bottom=0)
    axs[0].set_xlim(left=0)
    axs[0].axvline(p_c, ymax=f_max/plt.ylim()[1], ls = '--', c='red')
    axs[0].axhline(flow(p_c, p_max, V_0, n_green), xmax=p_c/plt.xlim()[1], ls = '--', c='red')
    axs[0].set_xlabel("Density in vehicle/km")
    axs[0].set_ylabel("Flux in vehicle/h")
    axs[0].set_title("Flux as a function of density - Greenshields model")
    axs[0].annotate(str("p_c = " + str(p_c)), (p_c,0), (p_c + 20, 150), arrowprops=dict(facecolor='black', shrink=0.01))
    axs[0].annotate(str("f_max = " + str(f_max)), (0,f_max), (20, f_max - 150), arrowprops=dict(facecolor='black', shrink=0.01))

    p = np.arange(0, p_max)
    v = V_0*(1-p/p_max)
    axs[1].plot(p,v)
    axs[1].set_xlabel("Density in vh/km")
    axs[1].set_ylabel("Speed in km/h")
    axs[1].set_title("Speed as a function of density - Greenshields model")

    plt.tight_layout()
    plt.show()

plot_greenshields(100, 50, 1)

