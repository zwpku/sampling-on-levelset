#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
from math import pi
import matplotlib.gridspec as gridspec

plt.figure(figsize=(8, 4))
axs = plt.gca()

N = int (raw_input('N='))

n_dof = N * (N-1) / 2

data_file = open('../working_dir/data/ex4_mcmc_v_norm_counter_%d.txt' % N, 'r')

n_conv, n_no_conv, size_s, v_b, n_bin = [ float(x) for x in data_file.readline().split() ]
bin_width = v_b / n_bin 

n_bin = int (n_bin)

counter_conv_data = [ float(x) for x in data_file.readline().split() ]

counter_no_conv_data = [ float(x) for x in data_file.readline().split() ]

tot_counter = [ counter_conv_data[i] + counter_no_conv_data[i] for i in range(n_bin) ]

counter_conv_data = [ x / (n_conv * bin_width) for x in counter_conv_data ]
#counter_no_conv_data = [ x / (n_no_conv * bin_width) for x in counter_no_conv_data ] 
tot_counter = [ x / ((n_conv + n_no_conv) * bin_width) for x in tot_counter ]

x_vec = np.linspace( bin_width * 0.5, v_b - bin_width*0.5, n_bin, endpoint=True )

axs.plot(x_vec, counter_conv_data , color='k',linewidth=1.5, linestyle =':', label=r'converge dist') 

#axs.plot(x_vec, counter_no_conv_data , color='b',linewidth=1.5, label=r'no converge dist') 

axs.plot(x_vec, tot_counter, color='r',linewidth=1.0, linestyle ='--', label=r'v dist') 

chi_density = [pow(x, n_dof-1) * exp(-(x/size_s) *(x/size_s)*0.5) for x in x_vec]

s = sum(chi_density) * bin_width
chi_density = [x / s for x in chi_density] 

axs.plot(x_vec, chi_density, color='k',linewidth=1.0, linestyle ='--', label=r'pdf Chi') 

axs.set_xlim(0, 5)
#axs.set_ylim(0.05, 0.35)
#axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
#axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
#axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
#axs.set_title(r'$\Theta$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

plt.show()

fig_file_name = '../fig/ex4_mcmc_v_dist_%d.eps' % N

savefig(fig_file_name, bbox_inches='tight')

