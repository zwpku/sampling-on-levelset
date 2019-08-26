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
mcmc_flag = int (raw_input('no-mcmc, mcmc, or both? (0/1/2)'))


if mcmc_flag == 0 or mcmc_flag == 2 :
    data_file = open('../working_dir/data/ex2_no_mcmc_counter_%d.txt' % N, 'r')
    n, xb, n_bin = [ int (x) for x in data_file.readline().split() ]
    bin_width = 2.0 * xb / n_bin 
    counter_data = [ (float(x) / (n * bin_width)) for x in data_file.readline().split() ]
    x_vec = np.linspace( -xb + bin_width * 0.5, xb - bin_width*0.5, n_bin, endpoint=True )
    axs.plot(x_vec, counter_data, color='k',linewidth=1.5, label=r'projection by $\Theta$') 

if mcmc_flag == 1 or mcmc_flag == 2 :
    data_file = open('../working_dir/data/ex2_mcmc_counter_%d.txt' % N, 'r')
    n, xb, n_bin = [ int (x) for x in data_file.readline().split() ]
    bin_width = 2.0 * xb / n_bin 
    counter_data = [ (float(x) / (n * bin_width)) for x in data_file.readline().split() ]
    x_vec = np.linspace( -xb + bin_width * 0.5, xb - bin_width*0.5, n_bin, endpoint=True )
    axs.plot(x_vec, counter_data, color='r',linestyle='--', marker='s',
            markevery=5, markersize=6, linewidth=1.5, label=r'MCMC') 

gaussian_density = [ exp(-x * x * 0.5) / sqrt(2 * pi) for x in x_vec]
axs.plot(x_vec, gaussian_density, color='k',linewidth=1.0, linestyle=':', label=r'std. Gaussian') 

axs.set_xlim(-xb, xb)
axs.set_ylim(-0.05, 0.5)
#axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
#axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
#axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
#axs.set_title(r'$\Theta$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(0.52, 0.68))

plt.show()

if mcmc_flag == 0:
    fig_file_name = '../fig/ex2_no_mcmc_trace_dist_%d.eps' % N
elif mcmc_flag == 1:
    fig_file_name = '../fig/ex2_mcmc_trace_dist_%d.eps' % N
else :
    fig_file_name = '../fig/ex2_both_trace_dist_%d.eps' % N

savefig(fig_file_name, bbox_inches='tight')
