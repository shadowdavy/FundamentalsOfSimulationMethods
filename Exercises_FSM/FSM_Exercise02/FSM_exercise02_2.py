import numpy as np
import matplotlib.pyplot as plt

k_B = 1.38 * 10**(-23)  # J / K
l_0 = 1 * 10**(-35)  # J m^3 s^-1
T_0 = 20000  # K
alpha = 10.0
beta = -0.5
n_H = 1 * 10**6  # m^⁻3
T_init = 1 * 10**7  # K


def rk2(f, y, t, dt):
    k1 = f(t, y)
    k2 = f(t + dt, y + k1 * dt)

    return (y + (1 / 2) * dt * (k1 + k2))

# y represents the temperature T
def func(dt, y):
    if (y <= T_0):
        l_T = l_0 * (y / T_0)**alpha
    else:
        l_T = l_0 * (y / T_0)**beta

    return (-(2 / (3 * k_B)) * n_H * l_T)
  
def loop(timestep):
    yrk = [T_init]
    T_break = T_init
    time = [0]
    steps = 0
    # loop until T = 6000 K
    while(T_break > 6000):
        yrk.append(rk2(func, yrk[-1], time[-1], timestep))
        time.append(time[-1] + timestep)
        T_break = yrk[-1]
        # print(T_break)
        steps += 1
    print(steps)
    yrk = np.array(yrk)
    time = np.array(time)
    return time, yrk

stepsize = 1 * 10**10  # s
data = loop(stepsize)
# plot
# plt.semilogy(time, yrk)
plt.semilogy(data[0], data[1], c = 'blue', label = "fixed timestep")
plt.title('temperature evolution T(t)')
plt.grid(True)

"""
data2 = loop(1 * 10**12)
plt.semilogy(data2[0], data2[1])
data2 = loop(1 * 10**11)
plt.semilogy(data2[0], data2[1])
data2 = loop(9.9 * 10**10)
plt.semilogy(data2[0], data2[1])
data2 = loop(9.7 * 10**10)
plt.semilogy(data2[0], data2[1])
"""

limit = 20  # K

def rk2_adabtive(f, y, t, dt):
    step1 = rk2(f, y, t, dt / 2)
    step2 = rk2(f, step1, t, dt / 2)
    full = rk2(f, y, t, dt)
    error = np.abs(full - step2)
    flag = True
    if (error > limit):
        dt = dt / 2
        # y = step2
        flag = False
    elif (10 * error < limit):
        dt = 2 * dt
        y = step2
    else:
        dt = dt
        y = full
    
    return y, dt, flag

yrk = [T_init]
T_break = T_init
time = [0]
steps = 0
stepsize = 1 * 10**10  # s
# loop until T = 6000 K
while(T_break > 6000):
    y_new, stepsize, flag = rk2_adabtive(func, yrk[-1], time[-1], stepsize)
    # print(stepsize)
    if (flag):
        yrk.append(y_new)
        time.append(time[-1] + stepsize)
        T_break = yrk[-1]
        #print(T_break)
    steps += 1
print(steps)
yrk = np.array(yrk)
time = np.array(time)
    
plt.semilogy(time, yrk, 'o', c = 'red', label = "adabtive timestep", ms = 2.5)
plt.semilogy(time, yrk, c = 'black', linewidth = 1.5)
# plt.title('temperature evolution T(t), adabtive dt')
# plt.grid(True)
plt.legend()

plt.savefig("FSM_exercise02_2.png")

plt.show()