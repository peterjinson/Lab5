import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import Orbit
#print(Orbit.orbit_traj(5,22,50,1,1)[0]) #orbit_traj(r0, v0, nStep, tau, NumericalMethod)
#thplot, rplot: orbital angle [rad] and radius [AU]

#print(np.shape(Orbit.orbit_traj(5,22,45,1,2)))

radius, theta = Orbit.orbit_traj(5,5,200,50,1)


posx = []
posy = []

fig, ax = plt.subplots()

artist = []

for t,r in zip(theta,radius):
    x_pos = r*np.cos(t)
    y_pos = r*np.sin(t)

    art = ax.scatter(x_pos,y_pos)

    artist.append(art)

ani = animation.ArtistAnimation(fig, artist, interval = 0.1, blit = True)
ani.save('Group4_Lab5_comet.m4a')


#for i in range(int(x)):

    #posx.append((Orbit.orbit_traj(5,2,50,1,2)[1][i])*np.cos(Orbit.orbit_traj(5,2,50,1,2)[0][i]))
    #posy.append((Orbit.orbit_traj(5,2,50,1,2)[1][i])*np.sin(Orbit.orbit_traj(5,2,50,1,2)[0][i]))


#print(len(posx))

#fig, ax = plt.subplots()
#ax.scatter(posx, posy)
#ax.grid(True)