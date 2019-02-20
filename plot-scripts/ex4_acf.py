#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
from math import pi
import matplotlib.gridspec as gridspec

plt.figure(figsize=(8, 6))
axs = plt.gca()

N = int (raw_input('N='))
mcmc_flag = int (raw_input('mcmc? (0/1)'))

if mcmc_flag == 1 :
    data_file = open('../working_dir/data/ex4_mcmc_acf_%d.txt' % N, 'r')
else :
    data_file = open('../working_dir/data/ex4_no_mcmc_acf_%d.txt' % N, 'r')

maxlag, mean_val, sigma, tau = [ float (x) for x in data_file.readline().split() ]
maxlag = int (maxlag)
  
acf_data = [ float(x) for x in data_file.readline().split() ]

axs.plot(acf_data, color='k',linewidth=1.5, label=r'acf') 

axs.set_xlim(0, 100)
axs.set_ylim(-0.1, 1.0)
#axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
#axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
axs.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
#axs.set_title(r'$\Theta$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

plt.show()

if mcmc_flag == 1:
    fig_file_name = '../fig/ex4_mcmc_acf_%d.eps' % N
else :
    fig_file_name = '../fig/ex4_no_mcmc_acf_%d.eps' % N

savefig(fig_file_name, bbox_inches='tight')
