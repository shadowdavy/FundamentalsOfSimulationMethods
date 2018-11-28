import numpy as np
import matplotlib.pyplot as plt

def rk2(f, y, t, dt):
    k1 = f(t, y)
    k2 = f(t + dt, y + k1 * dt)

    return (y + (1 / 2) * dt * (k1 + k2))

def rk4(f, y, t, dt):
    k1 = f(t, y)
    k2 = f(t + dt * 0.5, y + dt * 0.5 * k1)
    k3 = f(t + dt * 0.5, y + dt * 0.5 * k2)
    k4 = f(t + dt, y + dt * k3)

    return (y + (dt / 6) * (k1 + 2 * (k2 + k3) + k4))


M1 = 0.5
M2 = 1.0
L1 = 2.0
L2 = 1.0

def pendulum(t, y, m1 = M1, m2 = M2, l1 = L1, l2 = L2):
    
    phi1, phi2, q1, q2 = y

    g = 1
    f1 = (q1 - (l1 / l2) * np.cos(phi1 - phi2) * q2) / ((m1 + m2) * l1**2 - m2 * l2**2 * (np.cos(phi1 - phi2))**2)
    f2 = q2 / (m2 * l2**2) - (l1 / l2) * np.cos(phi1 - phi2) * ((q1 - np.cos(phi1 - phi2) * (q2 / l2)) / ((m1 + m2) * l1**2 - m2 * l1**2 * (np.cos(phi1 - phi2))**2))
    f3 = -m2 * l1 * l2 * f1 * f2 * np.sin(phi1 - phi2) - (m1 + m2) * g * l2 * np.sin(phi1)
    f4 = m2 * l1 * l2 * f1 * f2 * np.sin(phi1 - phi2) - m2 * g * l2 * np.sin(phi2)

    return np.array([f1, f2, f3, f4])

# initial conditions
phi1_0 = 50 * np.pi / 180 
phi2_0 = -120 * np.pi / 180 
phi1p_0 = phi2p_0 = 0  # phi point

y0 = np.array([phi1_0, phi2_0, phi1p_0, phi2p_0])


yrk4 = [y0]

dt = 0.05
time = np.arange(0, 100, dt)

# RK2
yrk2 = [y0]
for t in time:
    yrk2.append(rk2(pendulum, yrk2[-1], t, dt))
yrk2 = np.array(yrk2)

# RK4
yrk4 = [y0]
for t in time:
    yrk4.append(rk4(pendulum, yrk4[-1], t, dt))
yrk4 = np.array(yrk4)

fig = plt.figure()
ax = fig.add_subplot(211)
ax.scatter(yrk2[:,0], yrk2[:,2], s = 0.05)
ax.scatter(yrk2[:,1], yrk2[:,3], s = 0.05)
ax.set_title("RK 2")

ax2 = fig.add_subplot(212)
ax2.scatter(yrk4[:,0], yrk4[:,2], s = 0.05)
ax2.scatter(yrk4[:,1], yrk4[:,3], s = 0.05)
ax2.set_title("RK 4")

plt.savefig("FSM_exercise02_3.png")

plt.show()