from numpy import sqrt
import time
import RungeKuttaFehlberg as RKF
import numpy as np
# for Python2
# from tkinter import *   ## notice capitalized T in Tkinter
import scipy.integrate as integrate

import matplotlib.pyplot as plot
import matplotlib.animation as animation


class Orbit:
    GravConstant = 6.67408 * 10 ** (-11)
    M_e = 5.972 * 10 ** 24
    M_m =7.34767309*10**22
    h=0.1
    tol= 05e-14
    prevPositions = [[0],[384400000]]


    """

    Orbit Class

    init_state is [t0,x0,vx0,y0,vx0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """

    def __init__(self,
                 init_state,
                 G = GravConstant,
                 m1=M_e,
                 m2=M_m,
                 ):
        self.GravConst = G
        self.mPlanet1 = m1
        self.mPlanet2 = m2
        self.state = np.asarray(init_state, dtype='float')
        h1= self.h
        h2=self.h
        self.rkf54jorda = RKF.RungeKuttaFehlberg54(self.ydot, len(self.state[0]), h1, self.tol)
        self.rkf54månen = RKF.RungeKuttaFehlberg54(self.ydot, len(self.state[1]), h2, self.tol)
        self.prevPositions = self.prevPositions

    def getPos(self):
        return self.prevPositions

    def addPos(self, x,y):
        self.prevPositions[0].append(x)
        self.prevPositions[1].append(y)

    def position(self, i):
        """compute the current x,y positions of the pendulum arms"""
        x1 = self.state[i][1]
        y1 = self.state[i][3]
        return (x1, y1)

    # def energy(self):
    #     x = self.state[0][1]
    #     y = self.state[0][3]
    #     vx = self.state[0][2]
    #     vy = self.state[0][4]
    #     m1 = self.mPlanet1
    #     m2 = self.mSol
    #     G = self.GravConst
    #     U = -G * m1 * m2 / sqrt(x ** 2 + y ** 2)
    #     K = m1 * (vx ** 2 + vy ** 2) / 2
    #     return K + U

    def time_elapsed(self):
        return self.state[0][0], self.state[1][0]

    def step(self, i):
        w0 = self.state[i]
        if i ==0:
            self.state[i], E0 = self.rkf54jorda.safeStep(w0, i)
        else:
            self.state[i], E0 = self.rkf54månen.safeStep(w0, i)

    def ydot(self, x,i,h):
        if i == 0:
            m2 = self.mPlanet2
            px2 = self.state[1][1]
            py2 = self.state[1][3]

        elif i == 1:
            m2 = self.mPlanet1
            px2 = self.state[0][1]
            py2 = self.state[0][3]
        else :
            print("Feil i indeks fra state")
            exit(-1)

        px1 = x[1]
        py1 = x[3]
        vx1 = x[2]
        vy1 = x[4]
        z = np.zeros(5)
        dist = np.sqrt((px2-px1)**2 + (py2-py1)**2)
        z[0] = 1 #x[0] + h
        z[1] = vx1
        z[2] = (self.GravConst * m2 * (px2-px1))/(dist**3)
        # z[2] = self.GravConst * m2 / np.sqrt(px2 ** 2 + py2 ** 2) ** 3 *px2
        z[3] = vy1
        z[4] = (self.GravConst * m2 * (py2-py1))/(dist**3)
        # z[4] = self.GravConst * m2 / (np.sqrt(px2 ** 2 + py2 ** 2)) ** 3 * py2
        return z


# make an Orbit instance
#init_state is [t0,x0,vx0,y0,vx0],
orbit = Orbit([[0,0, 0, 0, 0],
                   [0.0,0, -1022, 384000000, 0]])
dt = 1. / 30 # 30 frames per second

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                       xlim=(-0.5*10**9, 0.5*10**9), ylim=(-0.5*10**9, 0.5*10**9))

trail, = axes.plot([], [], 'r--', lw=0.5)
lineA, = axes.plot([], [], 'o-b', lw=60, ms = 12)  # A blue planet 6*10**6
lineB, = axes.plot([], [], 'o-r', lw=17, ms = 3.4)  # A white planet

# line2, = axes.plot([], [], 'o-y', lw=2)  # A yellow sun
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
energy_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)


def init():
    """initialize animation"""
    lineA.set_data([], [])
    trail.set_data([], [])
    lineB.set_data([], [])
    # line2.set_data([], [])
    time_text.set_text('')
    energy_text.set_text('')
    return lineA, lineB, time_text, energy_text


def animate(i):
    """perform animation step"""
    global orbit, dt
    secondsPerFrame= 3600*24/30

    t0 = orbit.state[0][0]
    while (orbit.state[1][0] < t0 + secondsPerFrame):
        orbit.step(1)

    while (orbit.state[0][0] < t0 + secondsPerFrame):
        orbit.step(0)

    # if orbit.state[0][0] > 2629743 and orbit.state[1][0]>2629743 :
    pos = orbit.position(1)
    # print(orbit.state[0][0])

    x = pos[0]
    y = pos[1]
    orbit.addPos(x,y)
    trail.set_data(orbit.getPos())
    lineA.set_data(*orbit.position(0))
    lineB.set_data(*orbit.position(1))

    # line2.set_data([0.0, 0.0])
    t1, t2 = orbit.time_elapsed()
    antallDager = t1/(24*3600)

    time_text.set_text('time %.3f Days' % antallDager)
    # energy_text.set_text('energy = %.3f J' % orbit.energy())
    return lineA, lineB, time_text, energy_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 2000 * dt - (t1 - t0)

anim = animation.FuncAnimation(fig,  # figure to plot in
                               animate,  # function that is called on each frame
                               frames=900,  # total number of frames
                               interval=delay,  # time to wait between each frame.
                               repeat=False,
                               blit=True,
                               init_func=init  # initialization
                               )

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# plot.show()
