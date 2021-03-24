# Very blunt tool to test 2D polar-pattern simulation:

import numpy as np
import matplotlib.pyplot as plt

# define a number of points (N) to segment the circle around a microphone in multiple points 'n':
N = 249

# creating a circle:
nx = np.array([])
ny = np.array([])
nPoints = np.zeros((N, 2))

segments = np.linspace(0, 2, N)
for n in segments:
    x = np.cos(n*np.pi)
    if n < 1:
        y = np.sqrt(1 - (x*x))
    elif n >= 1:
        y = -np.sqrt(1 - (x*x))
    nx = np.append(nx, x)
    ny = np.append(ny, y)

# creating an array of 2D-vectors (dots around the circumference):
for i in range(0, N):
    vec = np.array([nx[i], ny[i]])
    nPoints[i] = vec

# prove it's a circle:
#plt.plot(*rays.T)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.show()

# copy the circle for every polar pattern:
cardioid = nPoints.copy()
biDir    = nPoints.copy()
omniDir  = nPoints.copy()
varDir   = nPoints.copy()

# scaling / shaping of the circles with |(1-s) + s(rayDirection * micDirection)|:
micDir = np.array([1, 0])
s_o = 0
s_c = 0.5
s_b = 1
s_v = 0.65
for n in range(0, N):
    omniDir[n]  *= abs(((1-s_o) + s_o * np.dot(omniDir[n], micDir)))
    cardioid[n] *= abs(((1-s_c) + s_c * np.dot(cardioid[n], micDir)))
    biDir[n]    *= abs(((1-s_b) + s_b * np.dot(biDir[n], micDir)))
    varDir[n]   *= abs(((1-s_v) + s_v * np.dot(varDir[n], micDir)))


# plots:
def render(signal, s):
    plt.plot(*signal.T, label="s = "+str(s))
    plt.legend(loc='best')
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.minorticks_on()
    plt.grid(True, linestyle='--')


fig = plt.figure()
ax = fig.gca()
ax.axhline(y=0, color='k', linewidth=.5)
ax.axvline(x=0, color='k', linewidth=.5)
ax.set_xticks(np.arange(-1, 1.1, 0.25))
ax.set_yticks(np.arange(-1, 1.1, 0.25))
plt.xlabel('x')
plt.ylabel('y')

render(omniDir, s_o)
render(cardioid, s_c)
render(biDir, s_b)
render(varDir, s_v)

plt.show()
