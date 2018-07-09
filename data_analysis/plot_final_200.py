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

filename1 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_1_3_8/Np_1_Nd_3_Hidden_3_cycle_*.dat'), key=lambda y: y.lower())
filename2 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_10_3_6/Np_10_Nd_3_Hidden_8_cycle_*.dat'), key=lambda y: y.lower())
filename3 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_50_3_6/Np_50_Nd_3_Hidden_6_cycle_*'), key=lambda y: y.lower())
filename4 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_70_3_6/Np_70_Nd_3_Hidden_6_cycle_*'), key=lambda y: y.lower())
filename5 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_90_3_6/Np_90_Nd_3_Hidden_6_cycle_*'), key=lambda y: y.lower())
filename6 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_100_3_3/Np_100_Nd_3_Hidden_3_cycle_*'), key=lambda y: y.lower())
filename7 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_120_3_3/Np_120_Nd_3_Hidden_3_cycle_*'), key=lambda y: y.lower())
filename8 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_150_3_3/Np_150_Nd_3_Hidden_3_cycle_*'), key=lambda y: y.lower())
filename9 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_180_3_3/Np_180_Nd_3_Hidden_3_cycle_*'), key=lambda y: y.lower())
filename10 = ns.natsorted(glob.glob('/home/dsacco/Desktop/thesis/Data_200_3_3/Np_200_Nd_3_Hidden_3_cycle_*'), key=lambda y: y.lower())

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

sigma=[];
z7=[]
Q7=[]
z7=np.asarray(z7)
z7.resize(1200)
Q7=np.asarray(Q7)
Q7.resize(1200)

sigma=[];
z8=[]
Q8=[]
z8=np.asarray(z8)
z8.resize(1200)
Q8=np.asarray(Q8)
Q8.resize(1200)

sigma=[];
z9=[]
Q9=[]
z9=np.asarray(z9)
z9.resize(1200)
Q9=np.asarray(Q9)
Q9.resize(1200)

sigma=[];
z10=[]
Q10=[]
z10=np.asarray(z10)
z10.resize(1200)
Q10=np.asarray(Q10)
Q10.resize(1200)


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

i=0;

for f in filename2:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
#    sigma[i]=block(x(axis=0));
#    y=block(x[:,1]);
#    z1[i]=y**.5
    y=stat.mean(x[:,1]);
    z2[i]=y
    u=block(x[:,1]);
    Q2[i]=u**.5;
    i=i+1;

i=0;

for f in filename3:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
#    sigma[i]=block(x(axis=0));
#    y=block(x[:,1]);
#    z1[i]=y**.5
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
    u=block(x[:,1]);
    Q5[i]=u**.5;
    y=stat.mean(x[:,1]);
    z5[i]=y
    i=i+1;

i=0;

for f in filename6:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    u=block(x[:,1]);
    Q6[i]=u**.5;
    y=stat.mean(x[:,1]);
    z6[i]=y
    i=i+1;

i=0;

for f in filename7:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    u=block(x[:,1]);
    Q7[i]=u**.5;
    y=stat.mean(x[:,1]);
    z7[i]=y
    i=i+1;

i=0;

for f in filename8:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    u=block(x[:,1]);
    Q8[i]=u**.5;
    y=stat.mean(x[:,1]);
    z8[i]=y
    i=i+1;

i=0;

for f in filename9:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    u=block(x[:,1]);
    Q9[i]=u**.5;
    y=stat.mean(x[:,1]);
    z9[i]=y
    i=i+1;

i=0;

for f in filename10:
    x=np.loadtxt(fname=f);
    x=x[len(x)-65536:len(x)]
    u=block(x[:,1]);
    Q10[i]=u**.5;
    y=stat.mean(x[:,1]);
    z10[i]=y
    i=i+1;

i=0;

#s=range(1200)
#s=np.asarray(s)
#plot1=plt.figure();
#plt.plot(s,z1,label='$\Delta t$=0.001')
##plt.plot(s,np.log10(z2),label='$\Delta t$=0.01')
##plt.plot(s,np.log10(z3),label='$\Delta t$=0.1')
##plt.plot(s,np.log10(z4),label='$\Delta t$=0.5')
##plt.plot(s,np.log10(z5),label='$\Delta t$=1.0')
##plt.plot(s,np.log10(z6),label='$\Delta t$=2.0')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final.pdf')
#pp.savefig(plot1)
#print(stat.mean(z1[400:len(z1)]))
#print(stat.mean(z1[350:len(z1)]))
#print(stat.mean(z1[300:len(z1)]))
#print(stat.mean(z1[250:len(z1)]))
#print(stat.mean(z1[200:len(z1)]))
#print(stat.mean(z1[150:len(z1)]))
#print(stat.mean(z1[100:len(z1)]))
#print(stat.mean(z1[50:len(z1)]))
#print("\007")
##pp.savefig(plot2)
##pp.savefig(plot3)
##pp.savefig(plot4)
##pp.savefig(plot5)
#pp.close()
#
#
#
#plot2=plt.figure();
#plt.plot(s,z2,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final2.pdf')
#pp.savefig(plot2)
#print(stat.mean(z2[400:len(z1)]))
#print(stat.mean(z2[350:len(z1)]))
#print(stat.mean(z2[300:len(z1)]))
#print(stat.mean(z2[250:len(z1)]))
#print(stat.mean(z2[200:len(z1)]))
#print(stat.mean(z2[150:len(z1)]))
#print(stat.mean(z2[100:len(z1)]))
#print(stat.mean(z2[50:len(z1)]))
#print("\007")
#pp.close()
#
#
#plot3=plt.figure();
#plt.plot(s,z3,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final3.pdf')
#pp.savefig(plot3)
#print(stat.mean(z3[400:len(z1)]))
#print(stat.mean(z3[350:len(z1)]))
#print(stat.mean(z3[300:len(z1)]))
#print(stat.mean(z3[250:len(z1)]))
#print(stat.mean(z3[200:len(z1)]))
#print(stat.mean(z3[150:len(z1)]))
#print(stat.mean(z3[100:len(z1)]))
#print(stat.mean(z3[50:len(z1)]))
#print("\007")
#pp.close()
#
#
#plot4=plt.figure();
#plt.plot(s,z4,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final4.pdf')
#pp.savefig(plot4)
#print(stat.mean(z4[400:len(z1)]))
#print(stat.mean(z4[350:len(z1)]))
#print(stat.mean(z4[300:len(z1)]))
#print(stat.mean(z4[250:len(z1)]))
#print(stat.mean(z4[200:len(z1)]))
#print(stat.mean(z4[150:len(z1)]))
#print(stat.mean(z4[100:len(z1)]))
#print(stat.mean(z4[50:len(z1)]))
#print("\007")
#pp.close()
#
#
##plot5=plt.figure();
##plt.plot(s,z5,label='$\Delta t$=0.001')
##plt.title('$\sigma$ evolution for various Importance sampling time step')
##plt.xlabel('SGD cycle')
##plt.ylabel('$\log10(\sigma)$')
###plt.ylim(-3.5, 0.5)
##plt.legend(loc=1);
##plt.show();
##pp = PdfPages('plot_final5.pdf')
##pp.savefig(plot5)
##print(stat.mean(z5[400:len(z1)]))
##print(stat.mean(z5[350:len(z1)]))
##print(stat.mean(z5[300:len(z1)]))
##print(stat.mean(z5[250:len(z1)]))
##print(stat.mean(z5[200:len(z1)]))
##print(stat.mean(z5[150:len(z1)]))
##print(stat.mean(z5[100:len(z1)]))
##print(stat.mean(z5[50:len(z1)]))
##print("\007")
##pp.close()
#
#
#plot6=plt.figure();
#plt.plot(s,z6,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final6.pdf')
#pp.savefig(plot6)
#print(stat.mean(z6[400:len(z1)]))
#print(stat.mean(z6[350:len(z1)]))
#print(stat.mean(z6[300:len(z1)]))
#print(stat.mean(z6[250:len(z1)]))
#print(stat.mean(z6[200:len(z1)]))
#print(stat.mean(z6[150:len(z1)]))
#print(stat.mean(z6[100:len(z1)]))
#print(stat.mean(z6[50:len(z1)]))
#print("\007")
#pp.close()
#
#
#plot7=plt.figure();
#plt.plot(s,z7,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final7.pdf')
#pp.savefig(plot7)
#print(stat.mean(z7[400:len(z1)]))
#print(stat.mean(z7[350:len(z1)]))
#print(stat.mean(z7[300:len(z1)]))
#print(stat.mean(z7[250:len(z1)]))
#print(stat.mean(z7[200:len(z1)]))
#print(stat.mean(z7[150:len(z1)]))
#print(stat.mean(z7[100:len(z1)]))
#print(stat.mean(z7[50:len(z1)]))
#print("\007")
#pp.close()
#
#
#plot8=plt.figure();
#plt.plot(s,z8,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final8.pdf')
#pp.savefig(plot8)
#print(stat.mean(z8[400:len(z1)]))
#print(stat.mean(z8[350:len(z1)]))
#print(stat.mean(z8[300:len(z1)]))
#print(stat.mean(z8[250:len(z1)]))
#print(stat.mean(z8[200:len(z1)]))
#print(stat.mean(z8[150:len(z1)]))
#print(stat.mean(z8[100:len(z1)]))
#print(stat.mean(z8[50:len(z1)]))
#print("\007")
#pp.close()
#
#
#plot9=plt.figure();
#plt.plot(s,z9,label='$\Delta t$=0.001')
#plt.title('$\sigma$ evolution for various Importance sampling time step')
#plt.xlabel('SGD cycle')
#plt.ylabel('$\log10(\sigma)$')
##plt.ylim(-3.5, 0.5)
#plt.legend(loc=1);
#plt.show();
#pp = PdfPages('plot_final9.pdf')
#pp.savefig(plot9)
#print(stat.mean(z9[400:len(z1)]))
#print(stat.mean(z9[350:len(z1)]))
#print(stat.mean(z9[300:len(z1)]))
#print(stat.mean(z9[250:len(z1)]))
#print(stat.mean(z9[200:len(z1)]))
#print(stat.mean(z9[150:len(z1)]))
#print(stat.mean(z9[100:len(z1)]))
#print(stat.mean(z9[50:len(z1)]))
#print("\007")
#pp.close()

print(stat.mean(z1[200:len(z1)]))
print(stat.mean(Q1[200:len(Q1)])/np.sqrt(len(Q1)-200))
print(stat.mean(z2[200:len(z2)]))
print(stat.mean(Q2[200:len(Q2)])/np.sqrt(len(Q2)-200))
print(stat.mean(z3[200:len(z3)]))
print(stat.mean(Q3[200:len(Q3)])/np.sqrt(len(Q3)-200))
print(stat.mean(z4[250:len(z4)]))
print(stat.mean(Q4[250:len(Q4)])/np.sqrt(len(Q4)-250))
print(stat.mean(z5[200:len(z5)]))
print(stat.mean(Q5[200:len(Q5)])/np.sqrt(len(Q5)-200))
print(stat.mean(z6[200:len(z6)]))
print(stat.mean(Q6[200:len(Q6)])/np.sqrt(len(Q6)-200))
print(stat.mean(z7[200:len(z7)]))
print(stat.mean(Q7[200:len(Q7)])/np.sqrt(len(Q7)-200))
print(stat.mean(z8[200:len(z8)]))
print(stat.mean(Q8[200:len(Q8)])/np.sqrt(len(Q8)-200))
print(stat.mean(z9[200:len(z9)]))
print(stat.mean(Q9[200:len(Q9)])/np.sqrt(len(Q9)-200))
print(stat.mean(z10[200:len(z10)]))
print(stat.mean(Q10[200:len(Q10)])/np.sqrt(len(Q10)-200))