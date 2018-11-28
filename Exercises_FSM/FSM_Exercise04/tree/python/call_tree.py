from tree import *
from random import *
import time
from copy import *
from math import sqrt
import matplotlib.pyplot as plt

#
# Create a set of randomly positioned particles
# For convenience we assume all masses to be 1.
# If we have in reality another mass, we can simply
# rescale our answers.
#
nparticles = 100
particles = []
for i in range(nparticles):
    x = random()
    y = random()
    z = random()
    particles.append([x,y,z])

#
# Now create the tree
#
q=TreeClass(particles)
q.insertallparticles()
q.computemultipoles(0)


print("starting tree gravity")
t0 = time.time()
q.allgforces(0.8)
t1 = time.time()
treegrav_dt = t1-t0
print("done in "+str(treegrav_dt)+" seconds\n")

fapprox = deepcopy(q.forces)

print("starting N^2 gravity")
t0 = time.time()

# ... TO BE FILLED IN ...
# print("Please replace this print statement with your own code")
#################
# appended code #
#################
q.forcesSum()

t1 = time.time()
fullgrav_dt = t1-t0
print("done in "+str(fullgrav_dt)+" seconds\n")

fexact = deepcopy(q.forces)
print("average particle node interaction: ", np.mean(q.nodecount))

# 
# Now compare the approximate and exact versions
#
# ... TO BE FILLED IN ...
#print("Please replace this print statement with your own code")

def absolut(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def error(f_tree, f_acc):
    return absolut(f_tree - f_acc) / absolut(f_acc)

colorArray = np.array(["red", "blue", "black"])

#nparticles = np.array([50, 100, 500, 1000])
nparticles = np.array([5, 10, 50])
angles = np.array([0.2, 0.4, 0.8])

fig = plt.figure()
meanPNI = np.zeros((len(angles), len(nparticles)))

for counter, angle in enumerate(angles):
    treegrav_dt = []
    fullgrav_dt = []
    eta = []
    etamean = []
    for j, npar in enumerate(nparticles):
        particles = []
        for i in range(npar):
            x = random()
            y = random()
            z = random()
            particles.append([x,y,z])

        q = TreeClass(particles)
        q.insertallparticles()
        q.computemultipoles(0)


        print("starting tree gravity")
        t0 = time.time()
        q.allgforces(angle)
        t1 = time.time()
        treegrav_dt.append(t1-t0)
        print("done in "+str(treegrav_dt[-1])+" seconds\n")
        fapprox = deepcopy(q.forces)

        print("starting N^2 gravity")
        t0 = time.time()
        q.forcesSum()
        t1 = time.time()
        fullgrav_dt.append(t1-t0)
        print("done in "+str(fullgrav_dt[-1])+" seconds\n")
        fexact = deepcopy(q.forces)
        
        qrange = len(q.forces)
        for i in range(qrange):
            eta.append(error(q.forces[i], q.accurateforces[i]))
        etamean.append(np.mean(eta))

        meanPNI[counter][j] = np.mean(q.nodecount)

    # print("etamean : ", etamean)
    print("average particle node interaction: ", meanPNI)
    #fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(nparticles, treegrav_dt, c = colorArray[counter], linewidth = 0.5)
    ax1.plot(nparticles, fullgrav_dt, c = colorArray[counter], linewidth = 0.5)
    ax1.scatter(nparticles, treegrav_dt, marker = '.', c = colorArray[counter], label = 'tree, angle: ' + str(angle))
    ax1.scatter(nparticles, fullgrav_dt, marker = 'x', c = colorArray[counter], label = 'direct')
    #ax1.set_xlabel('particles')
    ax1.set_ylabel('time in s')
    #plt.ylabel('time in s')
    plt.legend()

    ax2 = fig.add_subplot(212)
    ax2.plot(nparticles, etamean, c = colorArray[counter])
    ax2.scatter(nparticles, etamean, c = colorArray[counter], label = 'eta mean angle: ' + str(angle))
    ax2.set_ylabel('relative error')
    ax2.set_xlabel('particles')
    #plt.xlabel('particles')
    #plt.ylabel('relative error')
    #print(q.nodecount)
#plt.savefig("FSM_Exercise04.png")

plt.figure()
plt.plot(meanPNI)
plt.xlabel('angles')
plt.ylabel('particles')

plt.legend()
#plt.savefig("FSM_Exercise04_counter.png")
plt.show()
