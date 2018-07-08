##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Sun Jul  1 19:23:36 2018
#
#@author: dsacco
#"""
#
##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Sat Jun 30 17:59:05 2018
#
#@author: dsacco
#"""
#
#from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor, trim_zeros, append
#from numpy.linalg import inv
#from time import time
#from sys import argv
#
#import glob
#import numpy as np 
#import seaborn as sns
#import matplotlib.cm as cm
#import matplotlib.pyplot as plt
#import matplotlib.colors as colors
#from matplotlib.backends.backend_pdf import PdfPages
#import natsort as ns
#import statistics as stat
#
#def block(x):
#    # preliminaries
#    n = len(x); d = int(log2(n)); s, gamma = zeros(d), zeros(d);
#    mu = mean(x); t0 = time()
#
#    # estimate the auto-covariance and variances 
#    # for each blocking transformation
#    for i in arange(0,d):
#        n = len(x)
#        # estimate autocovariance of x
#        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
#        # estimate variance of x
#        s[i] = var(x)
#        # perform blocking transformation
#        x = 0.5*(x[0::2] + x[1::2])
#   
#    # generate the test observator M_k from the theorem
#    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]
#
#    # we need a list of magic numbers
#    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])
#
#    # use magic to determine when we should have stopped blocking
#    for k in arange(0,d):
#        if(M[k] < q[k]):
#            break
#    if(k >= d-1):
#        print("Warning: Use more data")
#    ans = s[k]/2**(d-k)
##    print("Runtime: %g sec" % (time()-t0))
##    print("Blocking Statistics :")
##    print("average            iterations      std. error")
##    print("%8g %20g %15g" % (mu, k, ans**.5))
#    return ans
#
## input data must be a power of two
##x = loadtxt(argv[1])
##x = x[int(0.75*len(x)):len(x)]
#
#sns.set();
#
#filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_1_3_8/Np_1_Nd_3_Hidden_3_cycle_*.dat'), key=lambda y: y.lower())
#
#
##filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_0.010000_cycle_*.dat'), key=lambda y: y.lower())
##filename3 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_0.100000_cycle_*.dat'), key=lambda y: y.lower())
##filename4 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_0.500000_cycle_*.dat'), key=lambda y: y.lower())
##filename5 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_1.000000_cycle_*.dat'), key=lambda y: y.lower())
##filename6 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Non_interacting/importance/Np_2_Nd_2_Hidden_2_Importance_timestep_2.000000_cycle_*.dat'), key=lambda y: y.lower())
#
#sigma=[];
#z1=[]
#Q1=[]
#z1=np.asarray(z1)
#z1.resize(1200)
#Q1=np.asarray(Q1)
#Q1.resize(1200)
#
#
#
#i=0;
#
#for f in filename1:
#    x=np.loadtxt(fname=f);
#    x=x[len(x)-65536:len(x)]
##    sigma[i]=block(x(axis=0));
##    y=block(x[:,1]);
##    z1[i]=y**.5
#    y=stat.mean(x[:,1]);
#    z1[i]=y
#    u=block(x[:,1]);
#    Q1[i]=y**.5;
#    i=i+1;
#
#
#print(stat.mean(z1[200:len(z1)]))
#print(stat.mean(Q1[200:len(Q1)])/np.sqrt(len(Q1)-200))
#
#s=range(1200)
#s=np.asarray(s)
#plot1=plt.figure();
#plt.plot(s,z1,label='Data')
##plt.plot(s,np.log10(z2),label='$\Delta t$=0.01')
##plt.plot(s,np.log10(z3),label='$\Delta t$=0.1')
##plt.plot(s,np.log10(z4),label='$\Delta t$=0.5')
##plt.plot(s,np.log10(z5),label='$\Delta t$=1.0')
##plt.plot(s,np.log10(z6),label='$\Delta t$=2.0')
#plt.title('Evolution of local energy for 1 particle in 3 dimensions')
#plt.xlabel('SGD cycle')
#plt.ylabel('Local energy')
#plt.ylim(2.0, 7.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final1_2.pdf')
#pp.savefig(plot1)
##pp.savefig(plot2)
##pp.savefig(plot3)
##pp.savefig(plot4)
##pp.savefig(plot5)
#pp.close()
#
#
#
#

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 19:23:36 2018

@author: dsacco
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_20_3_25/Np_20_Nd_3_Hidden_25_cycle_*.dat'), key=lambda y: y.lower())
#filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_30_3_27/Np_30_Nd_3_Hidden_27_cycle_*.dat'), key=lambda y: y.lower())

#filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_40_3_30/Np_40_Nd_3_Hidden_32_cycle_*.dat'), key=lambda y: y.lower())

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

#sigma=[];
#z2=[]
#Q2=[]
#z2=np.asarray(z2)
#z2.resize(1200)
#Q2=np.asarray(Q2)
#Q2.resize(1200)
#

i=0;

for f in filename1:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
#    sigma[i]=block(x(axis=0));
#    y=block(x[:,1]);
#    z1[i]=y**.5
    y=stat.mean(x[:,1]);
    z1[i]=y
    u=block(x[:,1]);
    Q1[i]=u**.5;
    i=i+1;

#i=0;
#
#for f in filename2:
#    x=np.loadtxt(fname=f);
#    x=x[len(x)-65536:len(x)]
##    sigma[i]=block(x(axis=0));
##    y=block(x[:,1]);
##    z1[i]=y**.5
#    y=stat.mean(x[:,1]);
#    z2[i]=y
#    u=block(x[:,1]);
#    Q2[i]=u**.5;
#    i=i+1;


print(stat.mean(z1[100:len(z1)]))
print(stat.mean(Q1[100:len(Q1)])/np.sqrt(len(Q1)-100))
print(stat.mean(z1[200:len(z1)]))
print(np.sqrt(np.sum(np.square(Q1[200:len(Q1)]))/(len(Q1)-200)**2))
print(stat.mean(z1[300:len(z1)]))
print(stat.mean(Q1[300:len(Q1)])/np.sqrt(len(Q1)-300))
print(stat.mean(z1[400:len(z1)]))
print(stat.mean(Q1[400:len(Q1)])/np.sqrt(len(Q1)-400))
print(stat.mean(z1[500:len(z1)]))
print(stat.mean(Q1[500:len(Q1)])/np.sqrt(len(Q1)-500))
print(stat.mean(z1[600:len(z1)]))
print(stat.mean(Q1[600:len(Q1)])/np.sqrt(len(Q1)-600))

#print(stat.mean(z2[100:len(z2)]))
#print(stat.mean(Q2[100:len(Q2)])/np.sqrt(len(Q2)-100))
#print(stat.mean(z2[200:len(z2)]))
#print(stat.mean(Q2[200:len(Q2)])/np.sqrt(len(Q2)-200))
#print(stat.mean(z2[300:len(z2)]))
#print(stat.mean(Q2[300:len(Q2)])/np.sqrt(len(Q2)-300))
#print(stat.mean(z2[400:len(z2)]))
#print(stat.mean(Q2[400:len(Q2)])/np.sqrt(len(Q2)-400))
#print(stat.mean(z2[500:len(z2)]))
#print(stat.mean(Q2[500:len(Q2)])/np.sqrt(len(Q2)-500))
#print(stat.mean(z2[600:len(z2)]))
#print(stat.mean(Q2[600:len(Q2)])/np.sqrt(len(Q2)-600))
#s=range(1200)
#s=np.asarray(s)
#plot1=plt.figure();
#plt.plot(s,z1,label='Data')
##plt.plot(s,np.log10(z2),label='$\Delta t$=0.01')
##plt.plot(s,np.log10(z3),label='$\Delta t$=0.1')
##plt.plot(s,np.log10(z4),label='$\Delta t$=0.5')
##plt.plot(s,np.log10(z5),label='$\Delta t$=1.0')
##plt.plot(s,np.log10(z6),label='$\Delta t$=2.0')
#plt.title('Evolution of local energy for 1 particle in 3 dimensions')
#plt.xlabel('SGD cycle')
#plt.ylabel('Local energy')
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_bonus1.pdf')
#pp.savefig(plot1)
##pp.savefig(plot2)
##pp.savefig(plot3)
##pp.savefig(plot4)
##pp.savefig(plot5)
#pp.close()
#
#
#plot2=plt.figure();
#plt.plot(s,z2,label='Data')
##plt.plot(s,np.log10(z2),label='$\Delta t$=0.01')
##plt.plot(s,np.log10(z3),label='$\Delta t$=0.1')
##plt.plot(s,np.log10(z4),label='$\Delta t$=0.5')
##plt.plot(s,np.log10(z5),label='$\Delta t$=1.0')
##plt.plot(s,np.log10(z6),label='$\Delta t$=2.0')
#plt.title('Evolution of local energy for 1 particle in 3 dimensions')
#plt.xlabel('SGD cycle')
#plt.ylabel('Local energy')
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_bonus2.pdf')
#pp.savefig(plot1)
##pp.savefig(plot2)
##pp.savefig(plot3)
##pp.savefig(plot4)
##pp.savefig(plot5)
#pp.close()

