from numpy import sqrt
import time
import RungeKuttaFehlberg as RKF
import numpy as np
# for Python2
# from tkinter import *   ## notice capitalized T in Tkinter
import scipy.integrate as integrate
import matplotlib.pyplot as plot
import matplotlib.animation as animation
import oppskytning
import saturn_v


class Orbit:
    GravConstant = 6.67408 * 10 ** (-11)
    M_e = 5.972 * 10 ** 24
    M_m = 7.34767309 * 10 ** 22
    h = 0.00000001
    tol = 05e-14
    prevPositions = [[0], [6371010]]

    """

    Orbit Class

    init_state is [t0,x0,vx0,y0,vx0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """

    def __init__(self,
                 init_state,
                 G=GravConstant,
                 m1=M_e,
                 m2=M_m,
                 ):
        self.GravConst = G
        self.mPlanet1 = m1
        self.state = np.asarray(init_state, dtype='float')
        self.rkf54 = RKF.RungeKuttaFehlberg54(self.ydot, len(self.state), self.h, self.tol)
        self.prevPositions = self.prevPositions

    def getPos(self):
        return self.prevPositions

    def addPos(self, x, y):
        self.prevPositions[0].append(x)
        self.prevPositions[1].append(y)

    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        x1 = self.state[1]
        y1 = self.state[3]
        x2 = self.state[5]
        y2 = self.state[7]
        return (x1, y1), (x2, y2)

    # def energy(self):
    #     pxJ = self.state[1]
    #     vxJ = self.state[2]
    #     pyJ = self.state[3]
    #     vyJ = self.state[4]
    #     pxM = self.state[5]
    #     vxM = self.state[6]
    #     pyM = self.state[7]
    #     vyM = self.state[8]
    #     mJorda = self.mPlanet1
    #     G = self.GravConst
    #     dist = np.sqrt((pxM - pxJ) ** 2 + (pyM - pyJ) ** 2)
    #     uTot = -G * mJorda * mManen / dist
    #     kJorda = mJorda * (vxJ ** 2 + vyJ ** 2) / 2
    #     kManen = mManen * (vxM **2 + vyM**2)/2
    #     return (kJorda + uTot + kManen )/(10**24)

    def time_elapsed(self):
        return self.state[0]

    def step(self):
        w0 = self.state
        self.state, E = self.rkf54.safeStep(w0)

    def ydot(self, x):
        mJorda = self.mPlanet1
        pxJ = x[1]
        vxJ = x[2]
        pyJ = x[3]
        vyJ = x[4]
        pxR = x[5]
        vxR = x[6]
        pyR = x[7]
        vyR = x[8]

        z = np.zeros(9)
        dist = np.sqrt((pxR - pxJ) ** 2 + (pyR - pyJ) ** 2)
        z[0] = 1
        z[1] = 0
        z[2] = 0
        z[3] = 0
        z[4] = 0
        z[5] = vxR
        z[6] = (-oppskytning.tyngdekraft((pyR-pyJ), saturn_v.masse(x[0])) - oppskytning.luftmotstand((pyR-pyJ), vyR, x[0]) + saturn_v.skyvekraft(x[0]))/oppskytning.masse(x[0])
        z[7] = vyR
        z[8] = (-oppskytning.tyngdekraft((pyR-pyJ), saturn_v.masse(x[0])) - oppskytning.luftmotstand((pyR-pyJ), vyR, x[0]) + saturn_v.skyvekraft(x[0]))/oppskytning.masse(x[0])
        # if(x[0]%1<0.01):print(x[0], oppskytning.tyngdekraft((pyR-pyJ), saturn_v.masse(x[0]))/saturn_v.masse(x[0]),oppskytning.luftmotstand((pyR-pyJ), vyR, x[0])/saturn_v.masse(x[0]), saturn_v.skyvekraft(x[0])/saturn_v.masse(x[0]))
        return z


# make an Orbit instance
# init_state: [t0, x0J, vx0J,  y0MJ   vy0J, x0R, vx0R,   y0R,    vy0R],
orbit = Orbit([0,   0,     0,     0,   0,    0,   0, 6371000, 0])
dt = 1. / 30  # 30 frames per second

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='auto', autoscale_on=False,
                       xlim=(-(2*6371000),2*6371000), ylim=(-(2*6371000),2*6371000))

trail, = axes.plot([], [], 'r--', lw=0.5)
lineA, = axes.plot([], [], 'o-b', lw=60, ms=128)  # A blue planet 6*10**6
lineB, = axes.plot([], [], 'o-r', lw=17, ms=3.4)  # A white planet

# line2, = axes.plot([], [], 'o-y', lw=2)  # A yellow sun
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
energy_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)


def init():
    """initialize animation"""
    lineA.set_data([], [])
    trail.set_data([], [])
    lineB.set_data([], [])
    time_text.set_text('')
    energy_text.set_text('')
    return lineA, lineB, time_text, energy_text


def animate(i):
    """perform animation step"""
    global orbit, dt
    secondsPerFrame = 4
    t0 = orbit.state[0]
    while orbit.state[0] < t0 + secondsPerFrame:
        orbit.step()
    # print("Fart: {}".format(orbit.state[8]))
    posJ, posR = orbit.position()
    # print(orbit.state)
    # print(orbit.rkf54.h)
    x = posR[0]
    y = posR[1]
    height = posR[1]-6371010
    orbit.addPos(x, y)
    trail.set_data(orbit.getPos())
    lineA.set_data(*posJ)
    lineB.set_data(*posR)
    t1 = orbit.time_elapsed()


    time_text.set_text('time %.3f S' % t1)
    energy_text.set_text('Height = %.5f m' % height)
    return lineA, lineB, time_text, energy_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 2000 * dt - (t1 - t0)

anim = animation.FuncAnimation(fig,  # figure to plot in
                               animate,  # function that is called on each frame
                               frames=300,  # total number of frames
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
# anim.save('saturn6sim.mp4', fps=30, extra_args=['-vcodec', 'libx264'])


plot.show()