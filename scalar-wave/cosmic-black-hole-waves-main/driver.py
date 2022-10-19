from helpers import *
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation

# The grid
def grid(metric):
    re = 1
    rc = 10
    L = 5

    if metric == 'dS':
        n = 251                  # spatial grid points
        x = np.linspace(0, L, n)     # the grid
        hgx = L/(n-1)      # spatial grid spacing
        hgt = 0.01               # time step

        endtime = 50            # end time
        tds = 2                 # downsampling in time
        nsteps = int(endtime/tds/hgt)  # number of time steps
        # Initial data
        x0 = 0.5*(L)
        sigma = (L)/40
        initial = np.zeros((2, len(x)))
        initial[1, :] = gaussian(x, x0, sigma)

    elif metric == 'SdS':  
        n = 251                  # spatial grid points
        x = np.linspace(re, rc, n)     # the grid
        hgx = (rc-re)/(n-1)      # spatial grid spacing
        hgt = 0.02               # time step

        endtime = 100            # end time
        tds = 2                 # downsampling in time
        nsteps = int(endtime/tds/hgt)  # number of time steps

        # Initial data
        x0 = 0.5*(re+rc)
        sigma = (re+rc)/40
        initial = np.zeros((2, len(x)))
        initial[1, :] = gaussian(x, x0, sigma)

    return initial, x,  hgx, hgt, tds, nsteps, L, re, rc




metric = 'dS'
coord = 'P'

initial, x,  hgx, hgt, tds, nsteps, L, re, rc = grid(metric)
print(L)
# Coefficients
ell = 0  # the spherical harmonic mode
coefs = set_coefs(x, ell, re, rc, L, coord, metric)

# Solve the wave equation
data = wave1dsolve(initial, coefs, hgx, hgt, tds, nsteps)

fig = plt.figure()

for frame_number in range(nsteps):
    plt.cla()
    plt.plot(x,data[frame_number,0,:])
    plt.pause(0.0005)
plt.show()


# def animate(frame_number):
#     plt.cla()
#     plt.plot(x,data[frame_number,0,:])
# anim = FuncAnimation(fig, animate, frames = 12*15)
# fig.suptitle(coord, fontsize=14)
# writervideo = animation.FFMpegWriter(fps=60) 
# anim.save(coord+'.mp4', writer=writervideo)
# plt.close()
