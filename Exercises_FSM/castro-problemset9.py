
# coding: utf-8

# In[133]:

get_ipython().magic('pylab inline')
import customvideo
from matplotlib import animation
get_ipython().magic("config InlineBackend.figure_formats=['svg']")
import numpy as np


# In[134]:

def riemann_roe_isothermal(rl, ul, vl, rr, ur, vr, cs):
    # assign values to U_R and U_L like the vector U is given on the sheet with U = (r, r*u, r*v)
    u_left = np.array([rl, rl * ul, rl * vl])
    u_right = np.array([rr, rr * ur, rr * vr])
    
    # do the same with vector F = (r*u, r*u^2 + r*c_s^2, r*u*s)
    fl = np.array([rl * ul, rl * ul**2 + rl * cs**2, rl * ul * vl])
    fr = np.array([rr * ur, rr * ur**2 + rr * cs**2, rr * ur * vr])
    
    # define ubar, vbar as given in the sheet
    ubar = (np.sqrt(rl) * ul + sqrt(rr) * ur) / (sqrt(rl) + sqrt(rr))
    vbar = (np.sqrt(rl) * vl + sqrt(rr) * vr) / (sqrt(rl) + sqrt(rr))
    
    # eigenvalues and eigenvector k
    eigenvalues = np.array([ubar - cs, ubar + cs, ubar])
    k = np.array([[1, ubar - cs, vbar], [1, ubar + cs, vbar], [0, 0, 1]])
    
    udiff = u_right - u_left 
    
    a = np.array([((ubar + cs) * udiff[0] - udiff[1]) / (2 * cs), 
                  -((ubar - cs) * udiff[0] + udiff[1]) / (2 * cs), 
                  udiff[2] - vbar * udiff[0]])

    # calculate the sum of vector F_star
    coefficients = a * np.abs(eigenvalues)
    k_coefficients = np.array([coefficients[0] * k[0], coefficients[1] * k[1], coefficients[2] * k[2]])
    
    # sum all components of eigenvector i elementwise (axis = 0), coefficients and eigenvalues included
    k_sum = np.sum(k_coefficients, axis = 0) 
    
    fs = []
    fs = 0.5 * (fl + fr) - 0.5 * k_sum
    
    return np.array(fs)


# In[135]:

cs = 2

# values stored as rl, ul, vl, rr, ur, vr
values = np.array([
        [1.0, 1.0, 2.0, 3.0, 1.0, 0.0],
        [2.5, 2.0, 3.0, 1.0, -3.0, -2.0],
        [2.0, -1.0, -2.0, 1.0, -1.0, 2.0]
    ])

for i in range(3):
    print(str(i + 1) + "th evaluation, ", riemann_roe_isothermal(values[i][0], values[i][1], values[i][2], values[i][3], values[i][4], values[i][5], cs))


# I could not verify the values given on the sheet in part d). Furthermore, I could not find the error which is responsible for the false values. I tried different variations of coding and array handling to check if the error occurs there. The last thing I could imagine is a wrong assignment of the vectors U_L and U_R, because I was not sure if I could just assign the given parameters of the function to this vector or if I need further calculation for getting U_L and U_R. Similarly for F_L, F_R.

# In[136]:

# this one need not to be evaluated, it is just a test of different array multiplications and sumations,
# because I thought the error is made there. But I guess the error occurs because of the wrong assignments of U_l and U_R.
"""
rl = 2.5
ul = 2
vl = 3
rr = 1
ur = -3
vr = -2
cs = 2

u_left = np.array([rl, rl * ul, rl * vl])
u_right = np.array([rr, rr * ur, rr * vr])
#print(u_right, u_left)    
fl = np.array([rl * ul, rl * ul**2 + rl * cs**2, rl * ul * vl])
fr = np.array([rr * ur, rr * ur**2 + rr * cs**2, rr * ur * vr])
#print(fl, fr)    
ubar = (np.sqrt(rl) * ul + sqrt(rr) * ur) / (sqrt(rl) + sqrt(rr))
vbar = (np.sqrt(rl) * vl + sqrt(rr) * vr) / (sqrt(rl) + sqrt(rr))
    
eigenvalues = np.array([ubar - cs, ubar + cs, ubar])
#print(eigenvalues)    
k = np.array([[1, ubar - cs, vbar], [1, ubar + cs, vbar], [0, 0, 1]])
udiff = u_right - u_left 

#print(k)
a = np.array([((ubar + cs) * udiff[0] - udiff[1]) / (2 * cs), -((ubar - cs) * udiff[0] + udiff[1]) / (2 * cs), udiff[2] - vbar * udiff[0]])

coefficients = a * np.abs(eigenvalues)
k_coefficients = np.array([coefficients[0] * k[0], coefficients[1] * k[1], coefficients[2] * k[2]])
#print(k_coefficients)
#print(np.sum(k_coefficients, axis = 0))
#sumk1 = coefficients[0] * np.array([1, ubar - cs, vbar])
#sumk2 = coefficients[1] * np.array([1, ubar + cs, vbar])
#sumk3 = coefficients[2] * np.array([0, 0, 1])
#print(sumk1)
#print(sumk2)
#print(sumk3)

#k[0] = coefficients[0] * k[0]
#k[1] = coefficients[1] * k[1]
#k[2] = coefficients[2] * k[2]

sumk = 0
for i in range(3):
    k[i] = coefficients[i] * k[i]
    sumk = sumk + k[i]
#print(sumk)
#sum_coefficients = np.array([udiff[0]*np.abs(eigenvalues[0]), udiff[1]*np.abs(eigenvalues[1]), udiff[2]*np.abs(eigenvalues[2])])
#sumall = sumk1 + sumk2 + sumk3
fs = 0.5 * (fl + fr) - 0.5 * sumk
print(fs)

#sum_coefficients = []
#for i in range(3):
#    sum_coefficients.append(np.sum([a, k[i]]))
#print(sum_coefficients)
"""


# In[ ]:




# In[ ]:



