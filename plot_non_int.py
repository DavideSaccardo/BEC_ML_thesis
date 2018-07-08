#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 16:24:11 2018

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

import natsort as ns
import statistics as stat

sns.set();
sns.set_style("whitegrid")
sns.set_context("paper")

filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/Np_2_Nd_2_Hidden_2_Metropolis_step_0.100000_cycle_*.dat'), key=lambda y: y.lower())
filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/Np_2_Nd_2_Hidden_2_Metropolis_step_0.010000_cycle_*.dat'), key=lambda y: y.lower())
filename3 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/Np_2_Nd_2_Hidden_2_Metropolis_step_0.500000_cycle_*.dat'), key=lambda y: y.lower())
filename4 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/Np_2_Nd_2_Hidden_2_Metropolis_step_1.000000_cycle_*.dat'), key=lambda y: y.lower())
filename5 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/Np_2_Nd_2_Hidden_2_Metropolis_step_2.000000_cycle_*.dat'), key=lambda y: y.lower())

sigma=[];
z1=[]
z1=np.asarray(z1)
z1.resize(101)
i=0;

z2=[]
z2=np.asarray(z2)
z2.resize(101)

z3=[]
z3=np.asarray(z3)
z3.resize(101)

z4=[]
z4=np.asarray(z4)
z4.resize(101)

z5=[]
z5=np.asarray(z5)
z5.resize(101)
#sns.set_context(rc={'lines.markeredgewidth': 1})

for f in filename2:
    x=np.loadtxt(fname=f);
    #print(np.size(x))
    x=x[len(x)-262144:len(x)]
    y=stat.mean(x[:,1]);
    z2[i]=y
    i=i+1;

i=0

for f in filename3:
    x=np.loadtxt(fname=f);
   # np.size(x)
    x=x[len(x)-262144:len(x)]
    y=stat.mean(x[:,1]);
    z3[i]=y
  #  print(i)
    i=i+1;

i=0

for f in filename4:
    x=np.loadtxt(fname=f);
    #np.size(x)
    x=x[len(x)-262144:len(x)]
    y=stat.mean(x[:,1]);
    z4[i]=y
    #print(i)
    i=i+1;

i=0

for f in filename5:
    x=np.loadtxt(fname=f);
    #np.size(x)
    x=x[len(x)-262144:len(x)]
    y=stat.mean(x[:,1]);
    z5[i]=y
    #print(i)
    i=i+1;

i=0

for f in filename1:
#    eval('loadtxt("/home/dsacco/Desktop/thesis/Non_interacting/Np_2_Nd_2_Hidden_2_Metropolis_step_0.010000_cycle_0.dat");');
#    eval('x_i=v_i[len(v_i)-262144:len(v_i)]');
#    eval('sigma_i=block(x_i(axis=1))');
#    eval('y_i=x_i.mean(axis=1)');
    x=np.loadtxt(fname=f);
    #np.size(x)
 #   print(x)
    x=x[len(x)-262144:len(x)]
   # print(len(x))
#    sigma[i]=block(x(axis=0));
    y=stat.mean(x[:,1]);
    z1[i]=y
   # print(z1[i]);
    #plt.plot(i,y[1])
    #print(i)
    i=i+1;

s=range(101)
s=np.asarray(s)
plot1=plt.figure();
plt.plot(s,z2,label='L=0.01')
#plot2=
plt.plot(s,z1,label='L=0.1')
#plot3=
plt.plot(s,z3,label='L=0.5')
#plot4=
plt.plot(s,z4,label='L=1.0')
#plot5=
plt.plot(s,z5,label='L=2.0')
plt.title('Local energy evolution for various BFM step lengths')
plt.xlabel('SGD cycle')
plt.ylabel('Local energy$\ [\hbar\omega_{ho}]$')
plt.legend();
plt.show();
pp = PdfPages('plot1.pdf')
pp.savefig(plot1)
#pp.savefig(plot2)
#pp.savefig(plot3)
#pp.savefig(plot4)
#pp.savefig(plot5)
pp.close()
#plt.xlim([-0.6, 0.6])
#plt.ylim(-0.1, 0.13)
#plt.plot(0,0.02)
    
#print(len(v))   
##v=v[0:512]
#x1=v[len(v)-1048576:len(v)]
##x1=v[len(v)-512:len(v)]
##x2=v[len(v)-262144:len(v)]
##x=v[len(v)-131072:len(v)]
#
##x1=x[len(x)-32768:len(x)]
##x2=x[len(x)-16384:len(x)]
##x=x[len(x)-65536:len(x)]
#
#print(len(x1))
##x = np.asarray(x)
##x2 = np.asarray(x2)
#x1 = np.asarray(x1)
#print(np.std(x1))
        

#plt.savefig("image.pdf")
#plt.savefig("image.png")
#plt.show()
