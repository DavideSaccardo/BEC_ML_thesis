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
from mpl_toolkits.axes_grid.inset_locator import inset_axes
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



sns.set();

#xdata=np.asarray([1,100,200,500,1000,  2000,5000,10000,15000,20000])
#ydata=np.asarray([2.414,2.66,2.86,3.30,3.84,4.61,6.12,7.76,8.98,9.98])
#
#y = func(xdata, 0.06586, 0.4848, 2.414)
#plot1=plt.figure();
#plt.scatter(xdata, ydata, label='Results')
#
#popt, pcov = curve_fit(func, xdata, ydata)
#
##plt.plot(xdata, func(xdata, *popt), 'r-', label='Fit') 
#         #: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
#f=interp1d(xdata,ydata, kind='nearest')
#plt.plot(xdata,f(xdata), 'r-', label='Interpolation')
#plt.title('Interpolation of the results')
#plt.xlabel('Number of Particles')
#plt.ylabel('$E_L/N_p$')
#plt.legend(loc=2);
#plt.show();
#pp = PdfPages('plot_fit.pdf')
#pp.savefig(plot1)
#pp.close()
#
#plot2=plt.figure();
##plt.plot(xdata, func(xdata, *popt), 'r-', label='Fit')
#f=interp1d(xdata,ydata,kind='nearest')
##f2 = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
#plt.plot(xdata,f(xdata), 'r-', label='Interpolation')
##plt.plot(xdata,f2(xdata), 'r-', label='Fit3')
#         #: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))


         
         
#fig, ax1 = plt.subplots()
#
## These are in unitless percentages of the figure size. (0,0 is bottom left)
#left, bottom, width, height = [0.6, 0.25, 0.3, 0.3]
#ax2 = fig.add_axes([left, bottom, width, height])
#
#

#24.443643344692994
#0.132869872786464
#49.31699614116059
#0.2629251903877939
#74.17151767252808
#0.3984364539634951
#99.15180194357605
#0.519307444132471
#124.76185252779084
#0.6467838484342862
#

N=np.asarray([10,20,30,40,50,80,100]);
E=np.asarray([2.44,2.47, 2.47,2.48,2.50,2.56,2.63]);
err=np.asarray([0.01,0.01,0.01, 0.01,0.01,0.01,0.02])
#ax1.errorbar(N, E, yerr=err,fmt='o',markersize=4, label='ML', color='blue')
NGP=np.array([100])
N_nonint=np.array([1])
EGP=np.array([2.66])
E_nonint=np.array([2.414])
#ax1.plot(NGP, EGP, label='GP', color='red')
#plt.xlim(0,110)
#plt.ylim(0,4)
#ax1.plot(N_nonint, E_nonint, label='Non interacting', color='green')
#plt.title('Comparison between ML and GP data')
#ax1.xlabel('Number of Particles')
#ax1.ylabel('$E_L/N_p$')
#ax1.legend(loc=2);
#ax2.plot(range(6)[::-1], color='green')
#plt.show();

sns.set_style("whitegrid")
sns.set_context("paper")
fig = plt.figure(); #figsize=(9, 4),facecolor='white')
ax = fig.add_subplot(111)
# the main axes is subplot(111) by default
plt.errorbar(N, E, yerr=err,fmt='o',markersize=2, label='ML', color='blue')
plt.plot(NGP, EGP, 'ro',markersize=2,label='GP')
plt.plot(N_nonint, E_nonint, 'g^',markersize=2,label='Non-interacting')
plt.xlim(0,110)
plt.ylim(0,4)
plt.title('Comparison between ML and GP data')
plt.xlabel('Number of Particles')
plt.ylabel('$E_L/N_p$$\ [\hbar\omega_{ho}]$')
plt.legend(loc=2);

# this is an inset axes over the main axes
inset_axes = inset_axes(ax,
                    width="100%", # width = 30% of parent_bbox
                    height=1.2, # height : 1 inch
                    bbox_to_anchor=(0.275,0.175,0.7,0.325), bbox_transform=ax.transAxes)

#plt.title('Probability')

#ax = fig.add_subplot(122)
# the main axes is subplot(111) by default
plt.errorbar(N, E, yerr=err,fmt='o',markersize=2, label='ML', color='blue')
NGP2=np.asarray([100,200,500])
EGP2=np.asarray([2.66,2.86,3.30])
plt.plot(NGP2, EGP2, 'ro',markersize=2,label='GP')
plt.plot(N_nonint, E_nonint, 'g^',markersize=2,label='Non interacting')
plt.axis([-1.5, 510, 0, 3.5])

pp = PdfPages('plot_finale.pdf')
pp.savefig(fig)
pp.close()