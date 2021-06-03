import numpy as np
import numpy.linalg as LA
import scipy.linalg
import matplotlib.pyplot as plt
import copy
import random
import warnings
import time
import math
import sys

warnings.simplefilter('ignore')


def lorenz96(x,F,J):
    """Lorenz 96 model."""
    # Compute state derivatives
    d = np.zeros(J)
    # First the 3 edge cases: i=1,2,N
    d[0] = (x[1] - x[J-2]) * x[J-1] - x[0]
    d[1] = (x[2] - x[J-1]) * x[0] - x[1]
    d[J-1] = (x[0] - x[J-3]) * x[J-2] - x[J-1]
    # Then the general case
    for i in range(2, J-1):
        d[i] = (x[i+1] - x[i-2]) * x[i-1] - x[i]
    # Add the forcing term
    d = d + F
    #print(F)
    # Return the state derivatives
    return d

def RK4(Var,dt,F,J):
    k1 = lorenz96(Var,F,J)
    k2 = lorenz96(Var+0.5*k1*dt,F,J)
    k3 = lorenz96(Var+0.5*k2*dt,F,J)
    k4 = lorenz96(Var+k3*dt,F,J)

    k = dt*(k1+2*k2+2*k3+k4)/6
    Var = Var + k
    return Var

def Mersenne(J):
    H = np.zeros(J)
    for i in range(J):
        a = random.normalvariate(0,1)
        H[i] = a
    return H

J = 40  # Number of variables
F = 8  # Forcing
epsilon = 0.001
F2 = 7.6
N = 7200
M = 40
L = 360
dt = 0.01
interval = 0.05 #unit
x = np.arange(40)

Var_t = np.loadtxt('Var_t.csv',delimiter=',',dtype='float') # true input
Var_o = np.loadtxt('Var_o.csv',delimiter=',',dtype='float') # obser input

b_list =  np.loadtxt('b_list.csv',delimiter=',',dtype='float')
Member = 24
ObNumM = 30
ObNumLast = 30
ObNum = 30
NumSt = 3
Inf = 1.5
DAWindow = 2 #日
T = math.floor(1460/DAWindow)
XxX = np.arange(8)

dataName =  str(Member)+"_"+str(Inf)+"_a_"+str(DAWindow)+".txt"
#dataName =  str(Member)+"_"+str(Inf)+"_"+str(DAWindow)+".txt"
f = open(dataName, mode='w')

Infs = np.array([1.8,2,2.5,3])
Infs = np.array([3.5,4,4.5,5])
i5M = np.size(Infs)

print("Member=",Member)
print("Obnum=",ObNumM)
print("T=",T*DAWindow)

Var_am = np.zeros((J,Member))
Var_a = copy.deepcopy(Var_t[0,:]) #+ Mersenne(J)
TT = 360
for i in range(TT*4):
    Var_a = RK4(Var_a,dt,F,J)
for i in range(Member):
    for ii in range(TT*4):
        Var_a = RK4(Var_a,dt,F,J)
    Var_am[:,i] = Var_a

print('ini finished')

I = np.eye(40)
aa = np.arange(-2,38)
aa2 = np.arange(40)
H = 1.0*I
R = 1.0*np.eye(NumSt)
Im = np.eye(Member)
delta = 10**(-3)

def GaspariCohn(r):
    if (r<0):
        r = -r
    if (np.abs(r)<=1):
        Gas = 1 - (1/4)*r**5 + (1/2)*r**4 + (5/8)*r**3 - (5/3)*r**2
    elif (np.abs(r)<=2):
        Gas = (1/12)*r**5 - (1/2)*r**4 + (5/8)*r**3 + (5/3)*r**2 - 5*r + 4 - (2/3)/r
    else:
        Gas = 0
    return Gas

def R_loc(R,NumSt):
    for i in range(NumSt):
        R[i,i] = R[i,i]/GaspariCohn(4*(i+1)/(NumSt+1)-2)
    return R

def LETKF(Member,delY_m_loc,R,delX_fm_loc,Im,xbar_f_loc,Var_o,H2):
    A1 = delY_m_loc@delY_m_loc.T + (Member-1)*R
    A2 = np.linalg.inv(A1)
    K = delX_fm_loc@delY_m_loc.T@A2 #(5)
    TTt = Im - delY_m_loc.T@A2@delY_m_loc

    D,U = LA.eig(TTt.T)
    Tx = U@(np.eye(Member)*scipy.linalg.sqrtm(np.eye(Member)*D))@U.T

    delX_am_loc = delX_fm_loc@Tx

    xbar_a_loc = xbar_f_loc + K@(Var_o - H2@xbar_f_loc)

    Var_am_loc =  delX_am_loc[1,:] + xbar_a_loc[1]

    return Var_am_loc

def Loc(H,location,b,J,aa):
    location2 = copy.deepcopy(location)
    for i in range(3):
        if (location[i]<0):
            location2[i] = location[i] + 40

    b2 = b - location[0]
    H4 = np.eye(3)
    H5 = np.delete(H4,b2,0)
    return H5

def LocR(R,location,b):
    b2 = b - location[0]
    R2 = np.delete(np.delete(R,b2,0),b2,1)
    return R2

Rini = R_loc(R,NumSt)

Var_amini = copy.deepcopy(Var_am)
for i5 in range(i5M):
    tt = np.arange(0,T*DAWindow,DAWindow)
    trPa = np.zeros(T)
    VarRMSE = np.zeros(T)
    Var_oRMSE = np.zeros(T)
    Var_aAvRMSE = np.zeros(J)
    b_l = np.zeros((4*DAWindow,J-ObNumM))
    dX_f_l = np.zeros((4*DAWindow,J,Member))
    Var_o_ii = np.zeros((4*DAWindow,J))
    M_j = np.zeros((J,J))
    Var_fm = np.zeros((J,Member))
    delX_fm = np.zeros((J,Member))
    delX_am = np.zeros((J,Member))
    Obseps = np.zeros((ObNumM,Member))
    Xbar_l = np.zeros((4*DAWindow,J))
    dX = np.zeros((J,Member))
    start = time.time()
    Inf = 1.4#copy.deepcopy(Infs[i5])
    Var_am = copy.deepcopy(Var_amini)
    for i in range(T):
        if(i%60==0):
            print(i*DAWindow,VarRMSE[i-1])
        Sum12 = np.zeros(Member)
        Sum13 = np.zeros((Member,Member))
        for i2 in range(4*DAWindow): #1day
            for i3 in range(5): #6hours
                for i4 in range(Member):
                    Var_am[:,i4] = RK4(Var_am[:,i4],dt,F,J)
            b_l[i2,:] = np.random.choice(aa,J-ObNumM,replace=False)
            #b_l[i2,:] = copy.deepcopy(b_list[4*i+i2,:])
            Var_o_ii[i2,:] = copy.deepcopy(Var_o[DAWindow*4*i+i2,:])
            Xbar_l[i2,:] = np.mean(Var_am,axis=1)
            for i3 in range(Member):
                dX_f_l[i2,:,i3] = Var_am[:,i3] - Xbar_l[i2,:] #dX

        Location = np.arange(-1,2)
        xbar_f = np.mean(Var_am,axis=1)
        for i2 in range(J):
            if (i2==39):
                Location = Location - 40
            for i3 in range(4*DAWindow):
                Var_o_loc = np.delete(Var_o_ii[i3,Location],b_l[i3,:]-Location[0])
                delX_fm_loc = dX_f_l[i3,Location,:]
                H6 = Loc(H,Location,b_l[i3,:],J,aa2)
                #R = scaling[i3]*Rini
                R = copy.deepcopy(Rini)
                Rloc = LocR(R,Location,b_l[i3,:])
                delY_m_loc = H6@delX_fm_loc
                aaaa = delY_m_loc.T@np.linalg.inv(Rloc)@delY_m_loc
                Sum13 = Sum13 + delY_m_loc.T@np.linalg.inv(Rloc)@delY_m_loc
                Sum12 = Sum12 + delY_m_loc.T@np.linalg.inv(Rloc)@(Var_o_loc - H6@Xbar_l[i3,Location])

            Pa_ny = np.linalg.inv((Member-1)*np.eye(Member)/Inf + Sum13) #(Member,Member)

            wbar = Pa_ny@Sum12 #(Member,1)

            D,U = LA.eig((Member-1)*Pa_ny)
            W_a = U@(np.eye(Member)*scipy.linalg.sqrtm(np.eye(Member)*D))@U.T
    #print(np.eye(Member)*scipy.linalg.sqrtm(np.eye(Member)*D))
            for i3 in range(Member):
                Var_am_sub = xbar_f[Location] + delX_fm_loc@(wbar + W_a[:,i3]) #(3,1)
                Var_am[i2,i3] = 1.0*Var_am_sub[1]

            Location = Location + 1

        meanVar_am = np.mean(Var_am,axis=1)
        for i2 in range(Member):
            dX[:,i2] = Var_am[:,i2]-meanVar_am
        Pa = (1/(Member-1))*dX@dX.T
        VarRMSE[i] = np.sqrt(((meanVar_am - Var_t[DAWindow*4*(i+1)-1,:])**2).mean())
        Var_oRMSE[i] = np.sqrt(((Var_o[DAWindow*4*(i+1)-1,:] - Var_t[DAWindow*4*(i+1)-1,:]) ** 2).mean())
        trPa[i] = np.sqrt(np.trace(Pa)/J)

    elapsed_time = time.time() - start

    Err2 = np.mean(VarRMSE[int(0.25*T):T])

    print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
    print("Err2=",Err2)
    print("Inf=",Inf)
    f.write("elapsed_time:{0}".format(elapsed_time) + "[sec]"+'\n')
    f.write("Err2="+str(Err2)+'\n')
    f.write("Inf="+str(Inf)+'\n')
    plt.figure(figsize=(28, 8), dpi=72)
    plt.plot(tt,trPa,label='trPa')
    plt.plot(tt,VarRMSE,label='analyticRMSE')
    #plt.ylim([0,7])
    plt.plot(tt,Var_oRMSE,label='observeRMSE')
    plt.legend()
    Name = str(Member)+"_"+str(Inf)+"_"+str(i5)+"_"+str(DAWindow)+".png"
    #Name = str(Member)+"_"+str(Inf)+"_"+str(i5)+"_"+str(DAWindow)+".png"
    plt.savefig(Name)
    #plt.show()
    plt.close()
f.close()
