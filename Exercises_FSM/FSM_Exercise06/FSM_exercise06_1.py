import numpy as np 
import matplotlib.pyplot as plt 
import copy

def tridiag(a, b, c, ka = -1, kb = 0, kc = 1):
    return np.diag(a, ka) + np.diag(b, kb) + np.diag(c, kc)

def matrixVectorMulti(matrix, vector):
    """
    result = []
    lenMatrix = len(matrix[0])
    lenVec = len(vector)
    if (lenMatrix != lenVec):
        print("length of matrix and vector are not equal")
        return -1
    else:
        for i in range(lenMatrix):
            total = 0
            for j in range(lenVec):
                total += matrix[i][j] * vec[j] 
            result.append(total)
        result = np.array(result)
        return result
    """
    return np.dot(matrix, vector) 

def febs(a, b, c, solutionVec):
    lenVec = len(solutionVec)
    a_copy = copy.deepcopy(a)
    b_copy = copy.deepcopy(b)
    c_copy = copy.deepcopy(c)
    vec = copy.deepcopy(solutionVec)

    # tridiagonal solver
    for i in range(1, lenVec):
        w = a_copy[i - 1] / b_copy[i - 1]
        b_copy[i] = b_copy[i] - w * c_copy[i - 1]
        vec[i] = vec[i] - w * vec[i - 1]

    # backward substitution
    x = b_copy
    x[-1] = vec[-1] / b_copy[-1]
    for i in range(lenVec - 2, -1, -1):
        x[i] = (vec[i] - c_copy[i] * x[i + 1]) / b_copy[i]

    return x


gridSize = 100
D = 0.5
epsilon = 1
T_0 = 1
L = 1
deltax = 2 * L / gridSize

diag_future = np.ones(gridSize - 1)
diag_present = -2 * np.ones(gridSize)
diag_past = np.ones(gridSize - 1)

# set matrix for boundary condition
diag_future[0] = 0
diag_past[-1] = 0
diag_present[0] = 1
diag_present[-1] = 1

A = tridiag(diag_past, diag_present, diag_future)

# text matrix multiplication
# print(diag_future, diag_present, diag_past)
# print(A)
# vec = np.array([2, 4, 5])
# print(matrixVectorMulti(A, vec))
# test elimination algorithm by comparing with numpy
# print(febs(diag_past, diag_present, diag_future, vec))
# print(np.linalg.solve(A, vec))

# apply to problem, epsilon is static
solutionVec = -(deltax**2 / D) * epsilon * np.ones(gridSize)
# boundary condition
solutionVec[0] = T_0
solutionVec[-1] = T_0  

result = febs(diag_past, diag_present, diag_future, solutionVec)
resultNp = np.linalg.solve(A, solutionVec)
# print(result)
# print(resultNp)
x = np.linspace(-L, L, gridSize)
plt.plot(x, result, color = 'red', label = 'grid size: ' + str(gridSize))
# plt.plot(x, resultNp)

# check if residual is almost zero
# print(solutionVec - matrixVectorMulti(A, result))
# why it is not exactly zero -> discretization and floating errors...

# change grid size
gridSize = 1000
deltax = 2 * L / gridSize
diag_future = np.ones(gridSize - 1)
diag_present = -2 * np.ones(gridSize)
diag_past = np.ones(gridSize - 1)

# set first future and last past value zero for fixed boundary condition
diag_future[0] = 0
diag_past[-1] = 0 
diag_present[0] = 1
diag_present[-1] = 1
A = tridiag(diag_past, diag_present, diag_future)
solutionVec = - deltax**2 / D * epsilon * np.ones(gridSize)
# boundary condition
solutionVec[0] = T_0
solutionVec[-1] = T_0  

result2 = febs(diag_past, diag_present, diag_future, solutionVec)
resultNp2 = np.linalg.solve(A, solutionVec)
# print(result2)
# print(resultNp2)
x = np.linspace(-L, L, gridSize)
plt.plot(x, result2, color = 'blue', label = 'grid size: ' + str(gridSize))
# plt.plot(x, resultNp2)


# Jacobi iteration
def jacobi(matrix, solutionVec, nit):
    lenVec = len(solutionVec)
    m_copy = copy.deepcopy(matrix)
    vec = copy.deepcopy(solutionVec)
    
    diag = np.diagflat(np.diag(m_copy))
    diagInv = np.linalg.inv(diag)
    diagUpper = -np.diagflat(np.diag(m_copy, 1), 1)
    diagLower = -np.diagflat(np.diag(m_copy, -1), -1)
    # x = np.zeros_like(vec)
    x = 2 * np.ones(len(vec))
    for i in range(nit):
        x = np.dot(diagInv, vec) + np.dot(diagInv, np.dot((diagLower + diagUpper), x))

    return x

gridSize = 8
deltax = 2 * L / gridSize
diag_future = np.ones(gridSize - 1)
diag_present = -2 * np.ones(gridSize)
diag_past = np.ones(gridSize - 1)

# set first future and last past value zero for fixed boundary condition
diag_future[0] = 0
diag_past[-1] = 0 
diag_present[0] = 1
diag_present[-1] = 1
A = tridiag(diag_past, diag_present, diag_future)
solutionVec = - deltax**2 / D * epsilon * np.ones(gridSize)
# boundary condition
solutionVec[0] = T_0
solutionVec[-1] = T_0  

n = 30
resultJac = jacobi(A, solutionVec, n)
# resultNpJac = np.linalg.solve(A, solutionVec)
# print(resultJac)
# print(resultNpJac)
x = np.linspace(-L, L, gridSize)
plt.plot(x, resultJac, color = 'black', label = 'Jacobi iteration, N: ' + str(gridSize))

# Jacobi N = 100
gridSize = 100
deltax = 2 * L / gridSize
diag_future = np.ones(gridSize - 1)
diag_present = -2 * np.ones(gridSize)
diag_past = np.ones(gridSize - 1)

# set first future and last past value zero for fixed boundary condition
diag_future[0] = 0
diag_past[-1] = 0 
diag_present[0] = 1
diag_present[-1] = 1
A = tridiag(diag_past, diag_present, diag_future)
solutionVec = - deltax**2 / D * epsilon * np.ones(gridSize)
# boundary condition
solutionVec[0] = T_0
solutionVec[-1] = T_0  

n = 30
resultJac = jacobi(A, solutionVec, n)
# resultNpJac = np.linalg.solve(A, solutionVec)
# print(resultJac)
# print(resultNpJac)
x = np.linspace(-L, L, gridSize)
plt.plot(x, resultJac, color = 'green', label = 'Jacobi iteration, N: ' + str(gridSize))


plt.xlabel('length [m]')
plt.ylabel('temperature [K]')
plt.legend()

plt.savefig('FSM_exercise06_1.png')

plt.show()  