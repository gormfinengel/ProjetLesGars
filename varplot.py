import numpy as np
import matplotlib.pyplot as plt

refttvar = [0., 0.00025, 0.0015, 0.0033, 0.019]
jumpttvar = [0., 0.00125, 0.0028, 0.0083, 0.0328]
dropttvar = [0., 0.00063, 0.0019, 0.004, 0.015]
sinrefttvar = [0., 0.00075, 0.0016, 0.0036, 0.0127]
sinjumpttvar = [0., 0.0009, 0.0023, 0.0075, 0.0229]
sindropttvar = [0., 0.00083, 0.0015, 0.0031, 0.0132]

var = np.array([refttvar, jumpttvar, dropttvar, sinrefttvar, sinjumpttvar, sindropttvar])

plt.plot(var[0], label="REFERENCE")
plt.plot(var[1], label="JUMP")
plt.plot(var[2], label="DROP")
plt.plot(var[3], label="SIN REFERENCE")
plt.plot(var[4], label="SIN JUMP")
plt.plot(var[5], label="SIN DROP")

plt.xticks([0., 1.,2.,3.,4.],[0.01, 0.15, 0.3, 0.5, 1.])
plt.xlabel("Noise variance (sigma^2)")
plt.ylabel("Travel time variance")
#plt.title("Travel times variances in relation to flux scenarios and noise variance")

plt.legend()
plt.show()