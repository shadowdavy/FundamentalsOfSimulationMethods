import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from hydro_iso_classic import *
from hydro_iso_riemann import *
from bc import *
#
# Model setup
#
bctype = 'periodic'
L      = 100.
csiso  = 1.
nx     = 100
nx = 10 * nx
cfl = 0.1
#cfl = 1
nt     = 400
t1     = 200.
dt     = cfl * t1/float(nt)
time   = np.arange(nt+1)*dt
x      = np.linspace(-L/2,L/2,nx)
dx     = x[1]-x[0]
x      = np.hstack((2*x[0]-x[1],x,2*x[-1]-x[-2]))
rho    = np.zeros_like(x) + 1.0 + np.exp(-x**2/200.)
rhov   = np.zeros_like(x)
cs     = np.zeros_like(x) + csiso

arr_rho  = np.zeros((nt+1,nx+2))
arr_rhov = np.zeros((nt+1,nx+2))
arr_rho[0,:]  = rho[:]
arr_rhov[0,:] = rhov[:]

#
# Time integration using first order scheme
#
q = np.vstack((rho,rhov))
implement_boundcond(q,bctype=bctype)
for it in range(nt):
    #q                = hydro_iso_classic_one_timestep(q,cs,dx,dt)
    q                = hydro_iso_riemann_one_timestep(q,cs,dx,dt)
    implement_boundcond(q,bctype=bctype)
    arr_rho[it+1,:]  = q[0,:]
    arr_rhov[it+1,:] = q[1,:]

#
# Animation of the result
#

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

count = 0
def update(frameNum, a0):
    global count,arr_rho,arr_rhov,x,nt
    y = arr_rho[count,:]
    a0.set_data(x, y)
    count = (count + 1) % nt
fig  = plt.figure()
ax   = plt.axes(xlim=(-50,50),ylim=(0, 2.4))
a0,  = ax.plot([], [])
anim = animation.FuncAnimation(fig, update, 
                               fargs=(a0,), 
                               interval=20)
anim.save('iso_standard_riemann_10.mp4', writer=writer)

plt.show()
