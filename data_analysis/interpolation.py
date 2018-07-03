#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 23:27:49 2018

@author: dsacco
"""

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor, trim_zeros, append
from numpy.linalg import inv
from time import time
from sys import argv
from scipy import interpolate

import glob
import numpy as np 
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import natsort as ns
import statistics as stat
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def func(x, a, b, c):
    return a*x**b + c

sns.set();

xdata=np.asarray([1,100,200,500,1000,  2000,5000,10000,15000,20000])
ydata=np.asarray([2.414,2.66,2.86,3.30,3.84,4.61,6.12,7.76,8.98,9.98])

y = func(xdata, 0.06586, 0.4848, 2.414)
plot1=plt.figure();
plt.scatter(xdata, ydata, label='Results')

popt, pcov = curve_fit(func, xdata, ydata)

#plt.plot(xdata, func(xdata, *popt), 'r-', label='Fit') 
         #: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
f=interp1d(xdata,ydata, kind='nearest')
plt.plot(xdata,f(xdata), 'r-', label='Interpolation')
plt.title('Interpolation of the results')
plt.xlabel('Number of Particles')
plt.ylabel('$E_L/N_p$')
plt.legend(loc=2);
plt.show();
pp = PdfPages('plot_fit.pdf')
pp.savefig(plot1)
pp.close()

plot2=plt.figure();
#plt.plot(xdata, func(xdata, *popt), 'r-', label='Fit')
f=interp1d(xdata,ydata,kind='nearest')
#f2 = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
plt.plot(xdata,f(xdata), 'r-', label='Interpolation')
#plt.plot(xdata,f2(xdata), 'r-', label='Fit3')
         #: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

#2.447648735836792
#0.04947335259741168
#24.443643344692994
#0.15634388833093144w
#49.31699614116059
#0.008285439253448332
#74.17151767252808
#0.01254700646707687

N=np.asarray([1,10,20,30,40,50]);
E=np.asarray([2.448,2.44,2.47, 2.47,2.48,2.50]);
err=np.asarray([0.01, 0.04,0.06,0.07, 0.08,0.09])
plt.errorbar(N, E, yerr=err,fmt='o',markersize=4, label='Data')
plt.title('Comparison between ML data and interpolation of results')
plt.xlabel('Number of Particles')
plt.ylabel('$E_L/N_p$')
plt.legend(loc=2);
plt.xlim(0,110)
plt.ylim(0,4)
plt.show();
pp = PdfPages('plot_fit2.pdf')
pp.savefig(plot2)
pp.close()