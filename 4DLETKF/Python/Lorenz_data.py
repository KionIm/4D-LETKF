import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import csv
import copy
import random
N = 40  # Number of variables
F = 8  # Forcing

def lorenz96(x):
    """Lorenz 96 model."""
    # Compute state derivatives
    d = np.zeros(N)
    # First the 3 edge cases: i=1,2,N
    d[0] = (x[1] - x[N-2]) * x[N-1] - x[0]
    d[1] = (x[2] - x[N-1]) * x[0] - x[1]
    d[N-1] = (x[0] - x[N-3]) * x[N-2] - x[N-1]
    # Then the general case
    for i in range(2, N-1):
        d[i] = (x[i+1] - x[i-2]) * x[i-1] - x[i]
    # Add the forcing term
    d = d + F
    #print(F)
    # Return the state derivatives
    return d

def RK4(Var,dt):
    k1 = lorenz96(Var)
    k2 = lorenz96(Var+0.5*k1*dt)
    k3 = lorenz96(Var+0.5*k2*dt)
    k4 = lorenz96(Var+k3*dt)

    k = dt*(k1+2*k2+2*k3+k4)/6
    Var = Var + k
    return Var


Var = F * np.ones(N)  # Initial state (equilibrium)
Var[20] += 0.008  # Add small perturbation to 20th variable
t = np.arange(0.0, 30.0, 0.01)
dt = 0.01
IM = int(2/0.05/8*8)
x = np.zeros(IM+1)
y = np.zeros(IM+1)
z = np.zeros(IM+1)
x[0] = Var[0]
y[0] = Var[1]
z[0] = Var[2]

def Pl(Var,i):
    Var2 = 5*(Var-8) + 0.25*(i+1)
    return Var2

fig = plt.figure()
x = np.arange(0,40)
Var2 = Pl(Var,-1)
plt.plot(x,Var2)

for i in range(365*20):
    Var = RK4(Var,dt)



def Mersenne(J):
    H = np.zeros(J)
    for i in range(J):
        a = random.normalvariate(0,1)
        H[i] = a
    return H

J = 40
epsilon_o = 1.0

ObNumM = 30
aa = np.arange(J)

f = open('Var_t.csv', 'w')
f2 = open('Var_o.csv', 'w')
f4 = open('b_list.csv', 'w')
writer = csv.writer(f, lineterminator='\n')
writer2 = csv.writer(f2, lineterminator='\n')
writer4 = csv.writer(f4, lineterminator='\n')

for ii in range(4*365*4):
    for i in range(5):
        Var = RK4(Var,dt)
    b = np.random.choice(aa,J-ObNumM,replace=False)
    writer.writerow(Var)
    ran_o = Mersenne(J)
    Var_o = Var + epsilon_o*ran_o
    writer2.writerow(Var_o)
    writer4.writerow(b)


f.close()
f2.close()
f4.close()

f3 = np.loadtxt('Var_t.csv',delimiter=',',dtype='float')

def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)


#adjustFigAspect(fig,aspect=10)
#ax = fig.add_subplot(111)

# Plot the first three variables
#fig = plt.figure()


#plt.xlabel('$x$')
#plt.ylabel('$time(days)$')
#plt.ylim([2.2,-0.2])
#plt.show()
