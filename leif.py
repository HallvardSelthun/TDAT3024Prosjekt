from numpy import sqrt
import time

import numpy as np
import scipy.integrate as integrate

import matplotlib.pyplot as plot
import matplotlib.animation as animation


class Orbit:
    """

    Orbit Class

    init_state is [t0,x0,vx0,y0,vx0],
    where (x0,y0) is the initial position
    , (vx0,vy0) is the initial velocity
    and t0 is the initial time
    """

    def __init__(self,
                 init_state,
                 G=1.0,
                 m1=1.0,
                 m2=1.0,
                 m3=1.0):
        self.GravConst = G
        self.mPlanet1 = m1
        self.mPlanet2 = m2
        self.mPlanet3 = m3
        self.state = np.asarray(init_state, dtype='float')
        # self.state2 = np.asarray(init_state2, dtype='float')
        # self.state3 = np.asarray(init_state3, dtype='float')


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
        return self.state[0][0]

    def step(self, h):
        """Uses the trapes method to calculate the new state after h seconds."""
        for i in range(3):
            x = self.state[i]
            s1 = self.ydot(x,i)
            s2 = self.ydot(x + h * s1,i)
            self.state[i] = x + h * (s1 + s2) / 2

    def ydot(self, x,i):
        G = self.GravConst
        if i == 0:
            m2 = self.mPlanet2
            px2 = self.state[1][1]
            py2 = self.state[1][3]
            m3 = self.mPlanet3
            px3 = self.state[2][1]
            py3 = self.state[2][3]
        elif i==1:
            m2 = self.mPlanet1
            px2 = self.state[0][1]
            py2 = self.state[0][3]
            m3 = self.mPlanet3
            px3= self.state[2][1]
            py3= self.state[2][3]
        else:
            m2 = self.mPlanet1
            px2 = self.state[0][1]
            py2 = self.state[0][3]
            m3 = self.mPlanet2
            px3 = self.state[1][1]
            py3 = self.state[1][3]
        px1 = x[1]
        py1 = x[3]
        vx1 = x[2]
        vy1 = x[4]
        dist1 = ((px2 - px1) ** 2 + (py2 - py1) ** 2)** (3/2)
        dist2 = ((px3 - px1) ** 2 + (py3 - py1) ** 2)** (3/2)
        z = np.zeros(5)
        z[0] = 1
        z[1] = vx1
        z[2] = (G*m2*(px2-px1)/dist1)+(G*m3*(px3-px1)/dist2)
        z[3] = vy1
        z[4] = (G*m2*(py2-py1)/dist1)+(G*m3*(py3-py1)/dist2)
        return z


# make an Orbit instance
#init_state is [t0,x0,vx0,y0,vx0],
orbit = Orbit([[0.0, -0.970, -0.466, 0.243, -0.433],
               [0.0, 0.970, -0.466, -0.243, -0.433],
               [0.0, 0.0, -2*(-0.466), 0.0,-2*(-0.433)]],
              1,1,1,1)
dt = 1. / 30  # 30 frames per second

# The figure is set
fig = plot.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                       xlim=(-2, 2), ylim=(-2, 2))

lineA, = axes.plot([], [], 'o-g', lw=2)  # A green planet
lineB, = axes.plot([], [], 'o-b', lw=2)  # A blue planetgreen
lineC, = axes.plot([], [], 'o-r', lw=2)  # A red planet

# line2, = axes.plot([], [], 'o-y', lw=2)  # A yellow sun
time_text = axes.text(0.02, 0.95, '', transform=axes.transAxes)
energy_text = axes.text(0.02, 0.90, '', transform=axes.transAxes)


def init():
    """initialize animation"""
    lineA.set_data([], [])
    lineB.set_data([], [])
    lineC.set_data([], [])
    # line2.set_data([], [])
    time_text.set_text('')
    energy_text.set_text('')
    return lineA,lineB, lineC, time_text, energy_text


def animate(i):
    """perform animation step"""
    global orbit, dt
    orbit.step(dt)
    lineA.set_data(*orbit.position(0))
    lineB.set_data(*orbit.position(1))
    lineC.set_data(*orbit.position(2))

    # line2.set_data([0.0, 0.0])
    time_text.set_text('time = %.1f' % orbit.time_elapsed())
    # energy_text.set_text('energy = %.3f J' % orbit.energy())
    return lineA, lineB,lineC, time_text, energy_text


# choose the interval based on dt and the time to animate one step
# Take the time for one call of the animate.
t0 = time.time()
animate(0)
t1 = time.time()

delay = 2000 * dt - (t1 - t0)

plot.show()