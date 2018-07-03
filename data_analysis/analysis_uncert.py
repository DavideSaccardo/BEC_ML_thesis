"""
Created on Sat Jun 30 17:59:05 2018

@author: dsacco
"""

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor, trim_zeros, append
from numpy.linalg import inv
from time import time
from sys import argv

import glob
import numpy as np 
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import natsort as ns
import statistics as stat

def block(x):
    # preliminaries
    n = len(x); d = int(log2(n)); s, gamma = zeros(d), zeros(d);
    mu = mean(x); t0 = time()

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if(k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
#    print("Runtime: %g sec" % (time()-t0))
#    print("Blocking Statistics :")
#    print("average            iterations      std. error")
#    print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans

# input data must be a power of two
#x = loadtxt(argv[1])
#x = x[int(0.75*len(x)):len(x)]

sns.set();

filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_1_3_3/Np_1_Nd_3_Hidden_3_cycle_*.dat'), key=lambda y: y.lower())
filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_10_3_8/Np_10_Nd_3_Hidden_8_cycle_*.dat'), key=lambda y: y.lower())
filename3 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_20_3_25/Np_20_Nd_3_Hidden_25_cycle_*.dat'), key=lambda y: y.lower())
filename4 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_30_3_27/Np_30_Nd_3_Hidden_27_cycle_*.dat'), key=lambda y: y.lower())
filename5 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_40_3_30/Np_40_Nd_3_Hidden_32_cycle_*.dat'), key=lambda y: y.lower())
filename6 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_50_3_40/Np_50_Nd_3_Hidden_40_cycle_*.dat'), key=lambda y: y.lower())

#filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_0.010000_cycle_*.dat'), key=lambda y: y.lower())
#filename3 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_0.100000_cycle_*.dat'), key=lambda y: y.lower())
#filename4 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_0.500000_cycle_*.dat'), key=lambda y: y.lower())
#filename5 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_1.000000_cycle_*.dat'), key=lambda y: y.lower())
#filename6 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_2.000000_cycle_*.dat'), key=lambda y: y.lower())

sigma=[];
z1=[]
Q1=[]
z1=np.asarray(z1)
z1.resize(1200)
Q1=np.asarray(Q1)
Q1.resize(1200)

sigma=[];
z2=[]
Q2=[]
z2=np.asarray(z2)
z2.resize(1200)
Q2=np.asarray(Q2)
Q2.resize(1200)

sigma=[];
z3=[]
Q3=[]
z3=np.asarray(z3)
z3.resize(1200)
Q3=np.asarray(Q3)
Q3.resize(1200)

sigma=[];
z4=[]
Q4=[]
z4=np.asarray(z4)
z4.resize(1200)
Q4=np.asarray(Q4)
Q4.resize(1200)

sigma=[];
z5=[]
Q5=[]
z5=np.asarray(z5)
z5.resize(1200)
Q5=np.asarray(Q5)
Q5.resize(1200)

sigma=[];
z6=[]
Q6=[]
z6=np.asarray(z6)
z6.resize(1200)
Q6=np.asarray(Q6)
Q6.resize(1200)


i=0;

for f in filename1:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    y=stat.mean(x[:,1]);
    z1[i]=y
    u=block(x[:,1]);
    Q1[i]=u**.5;
    i=i+1;

i=0;

for f in filename2:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    y=stat.mean(x[:,1]);
    z2[i]=y
    u=block(x[:,1]);
    Q2[i]=u**.5;
    i=i+1;

i=0;

for f in filename3:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    y=stat.mean(x[:,1]);
    z3[i]=y
    u=block(x[:,1]);
    Q3[i]=u**.5;
    i=i+1;
    
i=0;

for f in filename4:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    y=stat.mean(x[:,1]);
    z4[i]=y
    u=block(x[:,1]);
    Q4[i]=u**.5;
    i=i+1;

i=0;

for f in filename5:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    y=stat.mean(x[:,1]);
    z5[i]=y
    u=block(x[:,1]);
    Q5[i]=u**.5;
    i=i+1;
    

i=0;

for f in filename6:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    y=stat.mean(x[:,1]);
    z6[i]=y
    u=block(x[:,1]);
    Q6[i]=u**.5;
    i=i+1;
#
#print(stat.mean(z1[100:len(z1)]))
#print(stat.mean(Q1[100:len(Q1)])/np.sqrt(len(Q1)-100))

#Averages and uncertanties for N particles.
print(stat.mean(z1[200:len(z1)]))
print(np.sqrt(np.sum(np.square(Q1[200:len(Q1)]))/(len(Q1)-200)))
#print(block(Q1[len(Q1)-1024:len(Q1)])**.5)

print(stat.mean(z2[200:len(z2)]))
print(np.sqrt(np.sum(np.square(Q2[200:len(Q2)]))/(len(Q2)-200)))
#print(block(Q2[len(Q2)-1024:len(Q2)])**.5)

print(stat.mean(z3[200:len(z3)]))
print(np.sqrt(np.sum(np.square(Q3[200:len(Q3)]))/(len(Q3)-200)))
#print(block(Q3[len(Q3)-1024:len(Q3)])**.5)

print(stat.mean(z4[200:len(z4)]))
print(np.sqrt(np.sum(np.square(Q4[200:len(Q4)]))/(len(Q4)-200)))
#print(block(Q4[len(Q4)-1024:len(Q4)])**.5)

print(stat.mean(z5[200:len(z5)]))
print(np.sqrt(np.sum(np.square(Q5[200:len(Q5)]))/(len(Q5)-200)))
#print(block(Q5[len(Q5)-1024:len(Q5)])**.5)

print(stat.mean(z6[200:len(z6)]))
print(np.sqrt(np.sum(np.square(Q6[200:len(Q6)]))/(len(Q6)-200)))
#print(block(Q6[len(Q6)-1024:len(Q6)])**.5)

