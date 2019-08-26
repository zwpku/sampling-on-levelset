#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt

def xi(x1,x2):
    return 0.5 * (x1*x1/(c*c)+x2*x2-1) 

xb = 5
yb = 2
c = 3
scheme_id = 3 
data_file = open('../working_dir/data/ex1_traj_%d.txt' % (scheme_id) , 'r')
n = int (data_file.readline())
  
data = np.loadtxt(data_file, unpack=True)

y = np.array(data)

xi_vec = [ xi(y[0][i], y[1][i]) for i in range(n) ]

print "mean xi = %.4e,\tstd(xi) = %.4e\n" % (np.mean(xi_vec), np.std(xi_vec)) 

plt.figure(figsize=(8, 4))
ax = plt.gca()

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')

scatter(y[0], y[1], s=1, color='k')

plt.xlim(-xb, xb)
plt.ylim(-yb, yb)
xticks([-4, -2, 2, 4])
yticks([-2, 2])
ax.text(0.05, -0.32, '0', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)

fs = 30
#ax.set_title('$G_\epsilon(x)$', fontsize=fs)
ax.set_xlabel(r'$x_1$', fontsize=20)
ax.xaxis.set_label_coords(1.0, 0.49)
ax.set_ylabel(r'$x_2$', fontsize=20, rotation=0)
ax.yaxis.set_label_coords(0.525, 0.97)
ax.set_aspect('equal')

plt.tight_layout()
plt.show()
fig_file_name = '../fig/ex1_traj_%d.eps' % (scheme_id)
savefig(fig_file_name, bbox_inches='tight')
