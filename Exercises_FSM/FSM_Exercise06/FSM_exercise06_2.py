import numpy as np 
import matplotlib.pyplot as plt 

import numpy as np 
import matplotlib.pyplot as plt 
import copy

def tridiag(a, b, c, ka = -1, kb = 0, kc = 1):
    return np.diag(a, ka) + np.diag(b, kb) + np.diag(c, kc)
 

gridSize = 100
D = 0.5
epsilon = 1
T_0 = 1
L = 1
deltax = 2 * L / gridSize

diag_future = np.ones(gridSize - 1)
diag_present = -2 * np.ones(gridSize)
diag_past = np.ones(gridSize - 1)

# set first future and last past value zero for fixed boundary condition
diag_future[0] = 0
diag_past[-1] = 0 

A = tridiag(diag_past, diag_present, diag_future)