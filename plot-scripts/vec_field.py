#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
from math import pi
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

def xi(x,y):
    return 0.5*(x*x/(c*c) + y*y - 1)

lc = ['b', 'r', 'k', 'c', 'm', 'y']

fs = 30
# range of x, [-xb, xb] 
xb = 5
yb = 2
# number of discrete interval
n = 200
c = 3

fig, [ax1,ax2,ax3] = plt.subplots(1,3, figsize=(9,2), sharey=True)

Y, X = np.mgrid[-2:2:200j, -xb:xb:200j]
U = -X / (c * c) * xi(X,Y)
V = -Y * xi(X,Y)
speed = np.sqrt(U*U + V*V)

theta = np.linspace( 0, 2*pi, n, endpoint=True )
ax1.plot( c*np.cos(theta), np.sin(theta), color='k' )
lw = 5*speed / speed.max()
ax1.streamplot(X, Y, U, V, density=0.7, color='k', arrowsize=0.5, linewidth=0.8)
ax1.set_title(r'Flow map $\varphi$')
ax1.set_aspect('equal')
ax1.tick_params(length=0)
ax1.set_xticks(np.arange(-4, 5, step=2))
ax1.set_yticks(np.arange(-2, 3, step=2))

AA = 0.5
Y, X = np.mgrid[-2:2:200j, -xb:xb:200j]
U = -(X / (c * c) - AA * Y) * xi(X,Y)
V = -(AA * X / (c * c) + Y) * xi(X,Y)
speed = np.sqrt(U*U + V*V)

theta = np.linspace( 0, 2*pi, n, endpoint=True )
ax2.plot( c*np.cos(theta), np.sin(theta), color='k' )
lw = 5*speed / speed.max()
ax2.streamplot(X, Y, U, V, density=0.7, color='k', arrowsize=0.5, linewidth=0.8)
ax2.set_title(r'Flow map $\varphi^A$')
ax2.set_aspect('equal')
ax2.tick_params(length=0)
ax2.set_xticks(np.arange(-4, 5, step=2))
ax2.set_yticks(np.arange(-2, 3, step=2))

theta = np.linspace( 0, 2*pi, n, endpoint=True )
ax3.plot( c*np.cos(theta), np.sin(theta), color='k' )

line_len = 10
off = 0.5
theta_selected = np.linspace( 0, 2*pi, 35, endpoint=True )
thr_y = 0.1
for angle in theta_selected:
    pt = [c*np.cos(angle), np.sin(angle)]
    direction = [pt[0]/(c*c), pt[1] ]
    norm = np.sqrt(direction[0] * direction[0] + direction[1] * direction[1])
    direction = direction / norm
    pt1 = [pt[0] + line_len * direction[0], pt[1] + line_len * direction[1]]
    line=Line2D( [pt[0], pt1[0]], [pt[1], pt1[1]], color='k', linewidth=0.8)
    ax3.add_line(line)
    plt.arrow(pt[0] + off * (2-abs(pt[1])) * direction[0] , pt[1] + off *
            (2-abs(pt[1])) * direction[1], -0.2 * direction[0], -0.2 * direction[1], shape='full', lw=0.8,
            length_includes_head=True, head_width=.08)
    if abs(pt[1]) > thr_y :
        lam = - pt[1] / direction[1] 
        pt1 = [pt[0] + lam * direction[0], 0] 
        line=Line2D( [pt[0], pt1[0]], [pt[1], pt1[1]], color='k', linewidth=0.8)
        ax3.add_line(line)
        plt.arrow(pt[0] - off * direction[0] , pt[1] - off * direction[1], 0.1 * direction[0], 0.1 * direction[1], shape='full', lw=0.8, length_includes_head=True, head_width=0.08*abs(pt[1]))

plt.gcf().subplots_adjust(bottom=0.18) 
ax3.title.set_position([.5, 1.01])
ax3.set_title(r'Projection map $\Pi$')
ax3.set_aspect('equal')
ax3.set_xticks(np.arange(-4, 5, step=2))
ax3.tick_params(length=0)
ax3.set_yticks(np.arange(-2, 3, step=2))

plt.show()
plt.tight_layout()

fig_file_name = '../fig/vec_fields.eps' 
savefig(fig_file_name, bbox_inches='tight')

