import numpy as np
import matplotlib.pyplot as plt 

def scheme(u_old, a, deltax, deltat):
    u_new = u_old.copy()
    u_new[1:-1] = u_old[1:-1] - a * (deltat / (2 * deltax)) * (u_old[2:] - u_old[:-2])
    u_new[0] = u_old[0]
    u_new[-1] = u_new[-1]
    
    return u_new

def upwind(u_old, a, deltax, deltat):
    u_new = u_old.copy()
    u_new[1:-1] = u_old[1:-1] - a * (deltat / deltax) * (u_old[1:-1] - u_old[:-2])
    u_new[0] = u_old[0]
    u_new[-1] = u_new[-1]
    
    return u_new

def downwind(u_old, a, deltax, deltat):
    u_new = u_old.copy()
    u_new[1:-1] = u_old[1:-1] - a * (deltat / deltax) * (u_old[2:] - u_old[1:-1])
    u_new[0] = u_old[0]
    u_new[-1] = u_new[-1]
    
    return u_new

Nx = 100
L = 10
xini = np.linspace(0, L, Nx)
uini = np.zeros(Nx)
#uini[:Nx // 2] = 1.0
uini[np.where(xini < L / 2)] = 1

deltax = xini[1] - xini[0]
a = 1  # speed in advection equation
deltat = 0.03

def plotsequence(uini, timesteps, method = 'upwind'):
    u = uini.copy()
    for n in range(timesteps):
        if (method == 'scheme'):
            u = scheme(u, a, deltax, deltat)
        elif (method == 'upwind'):
            u = upwind(u, a, deltax, deltat)
        elif (method == 'downwind'):
            u = downwind(u, a, deltax, deltat)
        if n%5 == 0:
            plt.plot(xini, u, label = "t=%g"%(n * deltat))
    plt.legend(loc = "lower left")

timesteps = 100
fig = plt.figure()
plt.subplot(221)
plotsequence(uini, timesteps, 'scheme')
plt.title('standard scheme')
plt.subplot(222)
plotsequence(uini, timesteps, 'upwind')
plt.title('upwind')
plt.subplot(223)
plotsequence(uini, timesteps, 'downwind')
plt.title('downwind')

plt.savefig('FSM_exercise10_1_abc')

# change left boundary condition
uini[0] = 0.5

timesteps = 100
fig = plt.figure()
plt.subplot(221)
plotsequence(uini, timesteps, 'scheme')
plt.title('standard scheme, left boundary 0.5')
plt.subplot(222)
plotsequence(uini, timesteps, 'upwind')
plt.title('upwind, left boundary 0.5')
plt.subplot(223)
plotsequence(uini, timesteps, 'downwind')
plt.title('downwind, left boundary 0.5')

plt.savefig('FSM_exercise10_1_d')

# change right boundary condition
uini[0] = 1
uini[-1] = 0.5

timesteps = 100
fig = plt.figure()
plt.subplot(221)
plotsequence(uini, timesteps, 'scheme')
plt.title('standard scheme, right boundary 0.5')
plt.subplot(222)
plotsequence(uini, timesteps, 'upwind')
plt.title('upwind, right boundary 0.5')
plt.subplot(223)
plotsequence(uini, timesteps, 'downwind')
plt.title('downwind, right boundary 0.5')

plt.savefig('FSM_exercise10_1_e')

# 10 times smaller timesteps with 1000 iterations
uini[0] = 0.5
uini[-1] = 0.5
timesteps = 1000
deltax = 0.003
fig = plt.figure()
plt.subplot(221)
plotsequence(uini, timesteps, 'scheme')
plt.title('standard scheme, smaller timestep')
plt.subplot(222)
plotsequence(uini, timesteps, 'upwind')
plt.title('upwind, smaller timestep')
plt.subplot(223)
plotsequence(uini, timesteps, 'downwind')
plt.title('downwind, smaller timestep')

plt.savefig('FSM_exercise10_1_f')

# 10 times larger timesteps with 10 iterations
timesteps = 10
deltax = 0.3
fig = plt.figure()
plt.subplot(221)
plotsequence(uini, timesteps, 'scheme')
plt.title('standard scheme, larger timestep')
plt.subplot(222)
plotsequence(uini, timesteps, 'upwind')
plt.title('upwind, larger timestep')
plt.subplot(223)
plotsequence(uini, timesteps, 'downwind')
plt.title('downwind, larger timestep')

plt.savefig('FSM_exercise10_1_g')

plt.show()