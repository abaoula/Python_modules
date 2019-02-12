from __future__ import division,unicode_literals

from collections import Counter, OrderedDict
from scipy.integrate import odeint
from scipy.linalg import norm
from scipy.spatial.distance import euclidean
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import brewer2mpl
import random
import string
import os, sys
import warnings
from math import factorial
import pandas as pd
# import pytablewriter
from scipy import stats
import pandas as pd
from numba import jit

warnings.filterwarnings("ignore")


set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
set1 = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['lines.linewidth'] = '0.5'
plt.rcParams['axes.linewidth'] = '0.5'
plt.rcParams.update({'font.size': 8})

alphabet = list(string.ascii_lowercase)
alphabet1 = list(string.ascii_uppercase)

def cm2inch(value):
    return value/2.54


def binomial(x,p,n):
  le = len(x)
  p1 = np.empty((le))
  for i in range(le):
    p1[i] = factorial(n) / (factorial(x[i]) * factorial(n-x[i])) * p**x[i] * (1.-p)**(n-x[i])
  return [x, p1]

def table_to_markdown(def1):
  writer = pytablewriter.MarkdownTableWriter()
  writer.from_dataframe(df1)
  return writer.write_table()

def nonnan(a):
  return a[~np.isnan(a)]

def rsquare(x, y,y1):
  residuals = y - y1
  ss_res = np.sum(residuals**2)
  ss_tot = np.sum((y-np.mean(y))**2)
  r_squared = 1 - (ss_res / ss_tot)
  return r_squared

def histo(data, bins):
    x, y = np.histogram(data,np.linspace(data.min(), data.max(), bins+1))
    width= 0.95*(y[1]-y[0])
    return[x, y[:-1], 1.*x/np.sum(x), width]

def make_h(data,min, max, bins):
  p1 = np.zeros(bins)
  prob = np.linspace(min, max, bins+1)
  for i in range(len(data)):
    for j in range(bins):
      if data[i] >= prob[j] and data[i] < prob[j+1]:
        p1[j]= p1[j] + 1
  return p1

def count1(a):
  cc = Counter(a)
  c  = OrderedDict(sorted(cc.items()))
  return [np.asarray(c.keys()), np.asarray(c.values()), np.asarray(c.values())/(sum(np.asarray(c.values())))]


def hellinger(p, q):
    return norm(np.sqrt(p) - np.sqrt(q)) / _SQRT2

def pdist1(vec):
	l = len(vec)
	xv = np.unique(vec)
	l1 = len(xv)
	w1 = 1./l1
	count = np.zeros((l1,2))
	count[:,1] = xv
	for i in range(l):
		for j in range(l1):
			if np.round(vec[i], 5) == np.round(xv[j],5):
				count[j,0] = count[j,0]+ 1./l
	return [count, w1]

def pdistK(vec, n1):
	n = n1
	siz = len(vec)
	v =  np.linspace(-1, 1, n*(n-1)/2+1)
	l = len(v)
	w1 = 1./l
	count = np.zeros((l,2))
	count[:,1] = v
	for i in range(siz):
		aco =0
		for j in range(l):
			if np.round(vec[i],5) == np.round(v[j],5) and np.isnan(vec[i])==False:
				count[j,0] =count[j,0]+1./siz
	return [count, w1]

def finum(ax, x, y, text, color1):
  ax.text(x, y, text, transform=ax.transAxes, va='top', color=color1, )
  return ax


def figure(w1, ar, numbins):
	print 'width - height - # ticks'
	fig, ax = plt.subplots(ncols=1, nrows=1,figsize=(cm2inch(w1), cm2inch(ar)), dpi=100)
	ax.locator_params(nbins=numbins+1)
	#layout.cross_spines(ax=ax)

	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.margins(0.1, 0.1)
	#plt.tight_layout()

	return [fig, ax]

def hsubplot(w1, ar, col, numbins, nx, colo):
	print 'width - height - # cols - # ticks - sharey True or False'
	fig2, ax2 = plt.subplots(ncols=col, nrows=1,figsize=(cm2inch(w1), cm2inch(ar)), dpi=100, sharey=nx)
	plt.subplots_adjust(wspace=0.05, hspace=0.05)

	for i in range(0, col,1):
		# ax2[i].xaxis.set_ticks_position('bottom')
		# ax2[i].yaxis.set_ticks_position('left')
		ax2[i].locator_params(nbins=numbins+1)
		ax2[i].text(0.93, 0.96, alphabet[i], transform=ax2[i].transAxes, va='top', color=colo)
		ax2[i].set_xlabel('x')
		ax2[i].margins(0.1, 0.1)

	ax2[0].set_ylabel('y')
	#plt.tight_layout()

	return [fig2, ax2]

def vsubplot(w1, ar, row, numbins, nx):
	print 'width - height - # rows - # ticks - sharex True or False'
	fig3, ax3 = plt.subplots(ncols=1, nrows=row,figsize=(cm2inch(w1), cm2inch(ar)), dpi=100, sharex=nx)
	plt.subplots_adjust(wspace=0.05, hspace=0.05)
	for i in range(0, row,1):
		# ax3[i].xaxis.set_ticks_position('bottom')
		# ax3[i].yaxis.set_ticks_position('left')
		ax3[i].locator_params(nbins=numbins+1)
		#ax3[i].text(0.95, 0.95, alphabet[i], transform=ax3[i].transAxes, va='top', color='k')
		ax3[i].set_ylabel('y')
		ax3[i].margins(0.1, 0.1)

	ax3[row-1].set_xlabel('x')
	#plt.tight_layout()

	return [fig3, ax3]


def subplot(w1, ar, row, col, numbins, nx, ny, colo):
  print 'width - height - # row -  # col - # ticks - sharex True or False - sharey True or False'
  fig4, ax4 = plt.subplots(ncols=col, nrows=row,figsize=(cm2inch(w1), cm2inch(ar)), dpi=100, sharex=nx, sharey=ny)
  plt.subplots_adjust(wspace=0.06, hspace=0.06)
  ii = -1
  for i in range(0, row,1):
    for j in range(0,col,1):
      ii = ii+1
      # ax4[i,j].xaxis.set_ticks_position('bottom')
      # ax4[i,j].yaxis.set_ticks_position('left')
      ax4[i,j].margins(0.1, 0.1)
      ax4[i,j].locator_params(nbins=numbins+1)
      ax4[i,j].text(0.05, 0.96, alphabet[ii], transform=ax4[i,j].transAxes, va='top', color=colo)
  return [fig4, ax4]


print [key for key in locals().keys()
   if isinstance(locals()[key], type(sys)) and not key.startswith('__')]
