import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

def leapfrog(f, y, dt):
    #y[:6] = y[:6] + (dt / 2) * f(dt, y)[6:]
    #y[6:] = y[6:] + dt * y[:6]
    #y[:6] = y[:6] + (dt / 2) * f(dt, y)[6:]

    #y[6:] = y[6:] + (dt / 2) * f(dt, y)[:6]
    #y[:6] = y[:6] + dt * y[6:]
    #y[6:] = y[6:] + (dt / 2) * f(dt, y)[:6]

    #y[:6] = y[:6] + (dt / 2) * f(dt, y)[:6]
    #y[6:] = y[6:] + dt * y[:6]
    #y[:6] = y[:6] + (dt / 2) * f(dt, y)[:6]

    y_old = y
    x_half = y_old[6:] + (dt / 2) * y_old[:6]
    temp = np.concatenate((y_old[:6], x_half), axis = 0)
    v = y_old[:6] + dt * f(dt, temp)[:6]
    x = x_half + (dt / 2) * v
    y = np.concatenate((v, x), axis = 0)
   
    return y


def rk2(f, y, t, dt):
    k1 = f(t, y)
    k2 = f(t + dt, y + k1 * dt)

    return (y + (1 / 2) * dt * (k1 + k2))


def starplanet(t, y):
    G = 6.67384e-11
    ms = 1.99e30
    mp = 1e-3 * ms
    up, vp, wp, us, vs, ws, xp, yp, zp, xs, ys, zs = y
    g = G / ((xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2)**(3./2.)
    return np.array([- g * ms * (xp - xs),
                     - g * ms * (yp - ys),
                     - g * ms * (zp - zs),
                     - g * mp * (xs - xp),
                     - g * mp * (ys - yp),
                     - g * mp * (zs - zp),
                     up,
                     vp,
                     wp,
                     us,
                     vs,
                     ws])

# initial conditions
us0 = vs0 = ws0 = 0
xs0 = ys0 = zs0 = 0

xp0 = 1.496e11  # m
yp0 = zp0 = 0

up0 = 0
vp0 = 0.5 * 2.98e4  # m / s
wp0 = 0

y0 = np.array([up0, vp0, wp0, us0, vs0, ws0, xp0, yp0, zp0, xs0, ys0, zs0])

####################
# first few orbits #
####################

dt = 1.0 * 86400.
# five orbits
time = np.arange(0, 5 * 365 * 86400., dt)

# leapfrog
yleap = [y0]
for t in time:
    yleap.append(leapfrog(starplanet, yleap[-1], dt))
yleap = np.array(yleap)

fig = plt.figure()
ax = fig.add_subplot(221, projection = '3d')
ax.scatter(yleap[:,6], yleap[:,7], yleap[:,8], s = 0.05)
ax.scatter(yleap[:,9], yleap[:,10], yleap[:,11])
ax.set_title("LF 5 orbits")

# runge kutta 2
yrk = [y0]
for t in time:
    yrk.append(rk2(starplanet, yrk[-1], t, dt))
yrk = np.array(yrk)

ax3 = fig.add_subplot(223, projection = '3d')
ax3.scatter(yrk[:,6], yrk[:,7], yrk[:,8], s = 0.05)
ax3.scatter(yrk[:,9], yrk[:,10], yrk[:,11])
ax3.set_title("RK2 5 orbits")

##############
# 100 orbits #
##############
# 100 orbits
time = np.arange(0, 100 * 365 * 86400., dt)

# leapfrog
yleap = [y0]
for t in time:
    yleap.append(leapfrog(starplanet, yleap[-1], dt))
yleap = np.array(yleap)

ax2 = fig.add_subplot(222, projection = '3d')
ax2.scatter(yleap[:,6], yleap[:,7], yleap[:,8], s = 0.05)
ax2.scatter(yleap[:,9], yleap[:,10], yleap[:,11])
ax2.set_title("LF 100 orbits")

# runge kutta 2
yrk = [y0]
for t in time:
    yrk.append(rk2(starplanet, yrk[-1], t, dt))
yrk = np.array(yrk)

ax4 = fig.add_subplot(224, projection = '3d')
ax4.scatter(yrk[:,6], yrk[:,7], yrk[:,8], s = 0.05)
ax4.scatter(yrk[:,9], yrk[:,10], yrk[:,11])
ax4.set_title("RK2 100 orbits")

# save figure
# plt.savefig("FSM_exercise03_1.png")

plt.show()