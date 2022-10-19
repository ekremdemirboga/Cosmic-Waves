import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation

fig = plt.figure()
x,y = [], []

def animate(frame_number):
    plt.cla()
    x = np.linspace(0, 4, 1000)
    y = np.sin(2 * np.pi * (x - 0.01 * frame_number))
    plt.plot(x,y)

anim = FuncAnimation(fig, animate, interval=300)
fig.suptitle('Sine wave plot', fontsize=14)
writervideo = animation.FFMpegWriter(fps=60) 
anim.save('movie.mp4', writer=writervideo)
plt.close()