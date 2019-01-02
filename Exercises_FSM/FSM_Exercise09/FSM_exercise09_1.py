import numpy as np 
import matplotlib.pyplot as plt 


def randu(initialseed, iterations):
    if (initialseed % 2 == 0):
        print("initial value needs to be odd")
        return -1
    else:
        a = 65539
        m = 2147483648
        values = []
        values.append(np.int64(initialseed))

        for i in range(1, iterations):
            # compare to the first output values with wikipedia for I_0 = 0
            #if(i < 10):
            #    print(values[-1])
            values.append((a * values[-1]) % m)
            
        values = np.array(values)
        values = values / float(m)
        
        return values
"""
random_numbers = randu(1, 1000)
fig1 = plt.figure()
plt.scatter(random_numbers[::2], random_numbers[1::2])
plt.savefig("randu.png")
"""
# c++ file all numbers
fig1 = plt.figure()
file_r = np.loadtxt("randu.txt")
plt.scatter(file_r[::2], file_r[1::2], c = 'blue')

# c++ file restricted
fig2 = plt.figure()
file_rr = np.loadtxt("randuRestricted.txt")
plt.scatter(file_rr[:1000], file_rr[1000:], c = 'blue')
plt.xlim(0.2, 0.201)
plt.ylim(0.3, 0.301)
plt.savefig("randu_zoom.png")

# better random generator - restricted - draw individual points of uniform distribution
i = 0
ranx = []
rany = []
while (i < 500):
    ranx.append(np.random.uniform(0.2, 0.201))
    rany.append(np.random.uniform(0.3, 0.301))
    i = i + 1
ranx = np.array(ranx)
rany = np.array(rany)
fig3 = plt.figure()
plt.scatter(ranx, rany, c = 'red')
plt.xlim(0.2, 0.201)
plt.ylim(0.3, 0.301)
plt.savefig("random_zoom.png")

plt.show()







#########################################################
#########################################################
# not working correctly or long calculation time needed #
#########################################################

# long calculation time
"""
def randuRestricted(initialseed, iterations, x_interval = (0.2, 0.201), y_interval = (0.3, 0.301)):
    if (initialseed % 2 == 0):
        print("initial value needs to be odd")
        return -1
    else:
        a = 65539
        m = 2147483648
        values = []
        initial_2 = initialseed
        x_interval = np.array(x_interval) * m
        y_interval = np.array(y_interval) * m
        i = 0
        while (i < iterations):
            initial_1 = (a * initial_2) % m 
            initial_2 = (a * initial_1) % m

            if ((x_interval[0] < initial_1 < x_interval[1]) and (y_interval[0] < initial_2 < y_interval[1])): 
                values.append((initial_1 / m, initial_2 / m))
                i = i + 1
                print('number in intervall: ', i, "\n")
        return np.array(values)
"""

# not working
"""
points = 1000
while (True):
    random_numbers_zoom = randu(1, points)
    random_numbers_zoom = random_numbers_zoom[
                                            ((random_numbers_zoom[:] > 0.2) & (random_numbers_zoom[:] < 0.201)) |
                                            ((random_numbers_zoom[:] > 0.3) & (random_numbers_zoom[:] < 0.301))
                                            ]
    length = len(random_numbers_zoom)
    # print(length)
    if((length >= 900) and (length <= 1100)): break
    points = points * 2

rangex = random_numbers_zoom[random_numbers_zoom[:] < 0.201]
rangey = random_numbers_zoom[random_numbers_zoom[:] > 0.3]
if (len(rangex) < len(rangey)):
    rangey = rangey[:len(rangex)]
elif (len(rangey) < len(rangex)):
    rangex = rangex[:len(rangey)]

randomInteval = randuRestricted(1, 1000)
fig2 = plt.figure()
plt.scatter(*randomInteval.T)
#plt.scatter(rangex, rangey)
plt.xlim(0.1999, 0.2012)
plt.ylim(0.2999, 0.3012)
"""