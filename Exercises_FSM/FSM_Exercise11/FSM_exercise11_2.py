import numpy as np 
import matplotlib.pyplot as plt 

def advection(q_old, vel, deltax, deltat):
    q_new = q_old.copy()
    flux = np.zeros(len(vel))
 
    vpos = np.where(vel >= 0.)[0]
    vneg = np.where(vel < 0.)[0]
    flux[vpos] = q_old[vpos] * vel[vpos]
    flux[vneg] = q_old[vneg + 1] * vel[vneg]
    q_new[1:-1] = q_old[1:-1] - (deltat / deltax) * (flux[1:] - flux[:-1])
    q_new[0] = q_old[0]
    q_new[-1] = q_new[-1]

    return q_new

def momentumAdvection(momentum_old, r, cs, deltax, deltat):
    momentum_new = momentum_old.copy()
    momentum_new[1:-1] = momentum_old[1:-1] - cs**2 * (r[2:] - r[:-2]) / deltax

    return momentum_new

def hydro(rho_old, ve, cs, deltax, deltat):
    v = 0.5 * (ve[1:] + ve[:-1])
    rho = advection(rho_old, v, deltax, deltat)
    momentum = ve * rho
    momentum = advection(momentum, v, deltax, deltat)
    ve = momentumAdvection(momentum, rho, cs, deltax, deltat) / rho

    return rho, ve


Nx = 100
L = 50
xini = np.linspace(-L, L, Nx)
rho_0 = 1 + np.exp(- xini**2 / 200)
u_0 = np.zeros(Nx)
cs = 1

deltax = xini[1] - xini[0]
deltat = 0.5

t_start = 0
T = 50

fig1 = plt.figure()
while(t_start < T):
    rho_0, u_0 = hydro(rho_0, u_0, cs, deltax, deltat)

    if (round(0.5 * t_start) % 10 == 0):
        plt.plot(xini, rho_0)
    
    t_start += deltat
plt.savefig("FSM_exercise11_2b.png")

fig2 = plt.figure()
t_start = 0
while(t_start < T):
    rho_0, u_0 = hydro(rho_0, u_0, cs, deltax, deltat)

    if(round(0.5 * t_start) % 10 == 0):
        plt.plot(xini, rho_0)
    
    c = deltat * np.max(np.abs(u_0)) / deltax
    if (c >= 0.4):
        deltat = 0.4 * deltax / np.max(np.abs(u_0))
    
    t_start += deltat
plt.savefig("FSM_exercise11_2c.png")








plt.show()