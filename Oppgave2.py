import numpy as np
import RKFUnmodded as RKF
import math
import time
import matplotlib.pyplot as plot
from scipy.interpolate import spline

from datetime import datetime

"""
Velger oppgave 6.3.1d)

y1 = y1 + 3y2
y2 = 2y1 + 2y2
y1(0) = 5
y2(0) = 0
Med korrekt løsning
y1(t) = 3e**(−t) + 2*e**(4t)
y2(t) = −2*e**(−t) + 2*e**(4t)

"""
def fasit(t):
    y1 = 3*math.exp(-t) + 2*math.exp(4*t)
    y2 = -2*math.exp(-t) + 2*math.exp(4*t)
    return y1, y2

def F(y):
    z = np.ones(len(y))
    z[1]=y[1]+3*y[2]
    z[2]=2*y[1]+2*y[2]
    return z



h=0.1
accumulatedErrors =[]
globalErrors = [[],[]]
tolerances = []
times = []
tEnd=1.0
tol = (07e-16)/2

for i in range(40):
    W = np.array([0, 5, 0])
    tol =1.3*tol
    tolerances.append(tol)
    rkf54 = RKF.RungeKuttaFehlberg54(F, len(W), h, tol)
    accumulatedError = 0
    t0 = time.time()
    while (W[0]+rkf54.h < tEnd):
        W, E = rkf54.safeStep(W)
        accumulatedError +=E
    rkf54.setStepLength(tEnd - W[0])
    W, E = rkf54.step(W)
    t1 = time.time()
    timeElapsed = (t1-t0)
    times.append(1.0/timeElapsed)
    accumulatedErrors.append(accumulatedError)
    exact1, exact2 = fasit(tEnd)
    globalErrors[0].append(np.abs(exact1 - W[1]))
    globalErrors[1].append(np.abs(exact2 - W[2]))



for i in range(len(tolerances)):
    print("-----------------------")
    print(str(i+1)+'\t tol: '+str(tolerances[i]))
    print('Accumulated local error: '+str(accumulatedErrors[i]))
    print('Global error y1: '+str(globalErrors[0][i]))
    print('Global error y2: '+str(globalErrors[1][i]))
    print('Time elapsed calculating: '+str(times[i]))
    print("-----------------------")
fig = plot.figure()
# axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
#                        xlim=(0, tolerances[-1]), ylim=(0, times[-1]))

plot.loglog(tolerances, times)
plot.axis([tolerances[0], tolerances[-1],times[0], times[-1]])
plot.show()
