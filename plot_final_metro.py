#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 14:47:28 2018

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
    print("Runtime: %g sec" % (time()-t0))
    print("Blocking Statistics :")
    print("average            iterations      std. error")
    print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans

#f1="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_1_Nd_1_Hidden_1_Metropolis.dat";
#f2="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_1_Nd_2_Hidden_1_Metropolis.dat";
#f3="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_1_Nd_2_Hidden_2_Metropolis.dat";
#f4="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_1_Hidden_1_Metropolis.dat";
#f5="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_1_Hidden_2_Metropolis.dat";
#f6="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_1_Metropolis.dat";
#f7="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_2_Metropolis.dat";
#f8="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_3_Metropolis.dat";
#f9="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_4_Metropolis.dat";

#f1="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_1_Nd_1_Hidden_1_Importance.dat";
#f2="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_1_Nd_2_Hidden_1_Importance.dat";
#f3="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_1_Nd_2_Hidden_2_Importance.dat";
#f4="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_1_Hidden_1_Importance.dat";
#f5="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_1_Hidden_2_Importance.dat";
#f6="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_1_Importance.dat";
#f7="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_2_Importance.dat";
#f8="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_3_Importance.dat";
#f9="/home/dsacco/Desktop/thesis/Non_interacting/Final/Np_2_Nd_2_Hidden_4_Importance.dat";

#f1="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_1_Nd_3_Hidden_1_Metropolis.dat";
#f2="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_1_Nd_3_Hidden_2_Metropolis.dat";
#f3="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_1_Nd_3_Hidden_3_Metropolis.dat";
#f4="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_1_Metropolis.dat";
#f5="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_2_Metropolis.dat";
#f6="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_3_Metropolis.dat";
#f7="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_4_Metropolis.dat";
#f8="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_5_Metropolis.dat";
#f9="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_6_Metropolis.dat";

f1="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_1_Nd_3_Hidden_1_Importance.dat";
f2="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_1_Nd_3_Hidden_2_Importance.dat";
f3="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_1_Nd_3_Hidden_3_Importance.dat";
f4="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_1_Importance.dat";
f5="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_2_Importance.dat";
f6="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_3_Importance.dat";
f7="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_4_Importance.dat";
f8="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_5_Importance.dat";
f9="/home/dsacco/Desktop/thesis/Non_interacting/3rd/Np_2_Nd_3_Hidden_6_Importance.dat";



x=np.loadtxt(f1);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f2);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f3);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f4);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f5);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f6);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f7);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f8);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);

x=np.loadtxt(f9);
print(np.size(x))
x=x[len(x)-8388608:len(x)]
block(x[:,1]);
