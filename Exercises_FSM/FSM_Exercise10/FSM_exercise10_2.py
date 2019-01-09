import numpy as np 
import matplotlib.pyplot as plt 

def advection(q_old, vel, deltax, deltat):
    q_new = q_old.copy()
    # q_new[1:-1] = q_old[1:-1] - vel[:] * (deltat / deltax) * (q_old[1:-1] - q_old[:-2])
    lengthq = len(q_new)
    for i in range(1, lengthq - 1):
        # upwind
        if (vel[i] >= 0):
            q_new[i] = q_old[i] - (deltat / deltax) * (vel[i - 1] * q_old[i] - vel[i - 2] * q_old[i - 1])
        # downwind
        if (vel[i] < 0):
            q_new[i] = q_old[i] - (deltat / deltax) * (vel[i] * q_old[i + 1] - vel[i - 1] * q_old[i])
    q_new[0] = q_old[0]
    q_new[-1] = q_new[-1]

    return q_new

Nx = 100
L = 10
xini = np.linspace(0, L, Nx)
uini = np.zeros(Nx)
uini[np.where(xini < L / 2)] = 1

deltax = xini[1] - xini[0]
a = np.ones(Nx - 1)  # speed in advection equation
deltat = 0.03

timesteps = 100
fig = plt.figure()
u = uini.copy()
for n in range(timesteps):
    u = advection(u, a, deltax, deltat)
    if (n % 5 == 0):
        plt.plot(xini, u, label = "t=%g"%(n * deltat))
# plt.legend(loc = "lower left")
plt.savefig("FSM:exercise10_2_a.png")

# b)
v = -2 * (xini[1:-1] - (L / 2)) / L
v = np.pad(v, (1, 1), 'constant', constant_values = (v[-2], v[1]))
uini = np.zeros(Nx)
uini[np.where((np.abs(xini - (L / 2))) <= (L / 4))] = 1

timesteps = 100
fig = plt.figure()
u = uini.copy()
for n in range(timesteps):
    u = advection(u, v, deltax, deltat)
    if (n % 5 == 0):
        plt.plot(xini, u, label = "t=%g"%(n * deltat))
# plt.legend(loc = "lower left")
plt.savefig("FSM:exercise10_2_b.png")


plt.show()