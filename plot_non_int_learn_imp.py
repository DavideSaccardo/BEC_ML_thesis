#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 18:59:57 2018

@author: dsacco
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 23:46:55 2018

@author: dsacco
"""

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
import natsort as ns
import statistics as stat

sns.set();
sns.set_style("whitegrid")
sns.set_context("paper")

filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/learning/Np_2_Nd_2_Hidden_2_Importance_learn_0.050000_cycle_*.dat'), key=lambda y: y.lower())
filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/learning/Np_2_Nd_2_Hidden_2_Importance_learn_0.010000_cycle_*.dat'), key=lambda y: y.lower())
filename3 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/learning/Np_2_Nd_2_Hidden_2_Importance_learn_0.100000_cycle_*.dat'), key=lambda y: y.lower())
filename4 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/learning/Np_2_Nd_2_Hidden_2_Importance_learn_0.200000_cycle_*.dat'), key=lambda y: y.lower())
filename5 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/learning/Np_2_Nd_2_Hidden_2_Importance_learn_0.500000_cycle_*.dat'), key=lambda y: y.lower())
filename6 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/learning/Np_2_Nd_2_Hidden_2_Importance_learn_0.600000_cycle_*.dat'), key=lambda y: y.lower())

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

z6=[]
z6=np.asarray(z6)
z6.resize(101)

#sns.set_context(rc={'lines.markeredgewidth': 1})

for f in filename2:
    x=np.loadtxt(fname=f);
    #print(np.size(x))
    x=x[len(x)-262144:len(x)]
    y=stat.mean(x[:,1]);
    z2[i]=y
   # print(i)
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

for f in filename6:
    x=np.loadtxt(fname=f);
    #np.size(x)
    x=x[len(x)-262144:len(x)]
    y=stat.mean(x[:,1]);
    z6[i]=y
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
    y=x.mean(axis=0);
  #  print(len(y));
   # print(y[1])
    z1[i]=y[1]
   # print(z1[i]);
    #plt.plot(i,y[1])
    #print(i)
    i=i+1;

s=range(101)
s=np.asarray(s)
plot1=plt.figure();
plt.plot(s,z2,label='$\eta$=0.01')
#plot2=
plt.plot(s,z1,label='$\eta$=0.05')
#plot3=
plt.plot(s,z3,label='$\eta$=0.10')
#plot4=
plt.plot(s,z4,label='$\eta$=0.20')
#plot5=
plt.plot(s,z5,label='$\eta$=0.50')
plt.plot(s,z6,label='$\eta$=0.60')
plt.title('Local energy evolution for various ISMH learning rates')
plt.xlabel('SGD cycle')
plt.ylabel('Local energy$\ [\hbar\omega_{ho}]$')
plt.legend();
plt.show();
pp = PdfPages('plot6.pdf')
pp.savefig(plot1)
#pp.savefig(plot2)
#pp.savefig(plot3)
#pp.savefig(plot4)
#pp.savefig(plot5)
pp.close()