import matplotlib.pyplot as plt
import numpy as np

def flow(p, p_max, V_0, n_green):
    return V_0*p*(1-(p*(1/p_max)))

def plot_greenshields(p_max, V_0, n_green):
    p = np.arange(0, p_max)
    p_c = (p_max)/(n_green+1)**(1/n_green)
    f = flow(p, p_max, V_0, n_green)
    f_max = flow(p_c, p_max, V_0, n_green)
    plt.plot(p,f)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.axvline(p_c, ymax=f_max/plt.ylim()[1], ls = '--', c='red')
    plt.axhline(flow(p_c, p_max, V_0, n_green), xmax=p_c/plt.xlim()[1], ls = '--', c='red')
    plt.xlabel("Densité en vh/km")
    plt.ylabel("Flux en vh/h")
    plt.title("Distribution du flux en fonction de la densité - Modèle de Greenshields")
    plt.annotate(str("p_c = " + str(p_c)), (p_c,0), (p_c + 20, 150), arrowprops=dict(facecolor='black', shrink=0.01))
    plt.annotate(str("f_max = " + str(f_max)), (0,f_max), (20, f_max - 150), arrowprops=dict(facecolor='black', shrink=0.01))
    plt.show()

plot_greenshields(150, 50, 1)

