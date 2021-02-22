import numpy as np
import matplotlib.pyplot as plt

rayAmount = 249
sigX = np.array([])
sigY = np.array([])

# creating a circle:
k = np.linspace(0, 2, rayAmount)
for n in k:
    x = np.cos(n*2*np.pi)
    if n < 1:
        y = np.sqrt(1 - (x*x))
    elif n >= 1:
        y = -np.sqrt(1 - (x*x))
    sigX = np.append(sigX, x)
    sigY = np.append(sigY, y)

# prove it's a circle:
# plt.plot(sigX[0:rayAmount], sigY[0:rayAmount])
# plt.gca().set_aspect('equal', adjustable='box')
# plt.xlabel("x")
# plt.ylabel("y")

# creating an array of 2D-vectors (dots around the circumference):
rays = np.zeros((rayAmount, 2))
for i in range(0, rayAmount):
    vec = np.array([sigX[i], sigY[i]])
    rays[i] = vec

# test, whether those dots are still a circle:
#plt.plot(*rays.T)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.show()

# copy the circle for every polar pattern:
cardioid = rays.copy()
biDir = rays.copy()
omniDir = rays.copy()
varDir = rays.copy()

# scaling / shaping of the circles with |(1-s) + s(rayDirection * micDirection)|:
micDir = np.array([1, 0])
s_o = 0
s_c = 0.5
s_b = 1
s_v = 0.65
for n in range(0, rayAmount):
    omniDir[n] *= abs(((1-s_o) + s_o * np.dot(omniDir[n], micDir)))
    cardioid[n] *= abs(((1-s_c) + s_c * np.dot(cardioid[n], micDir)))
    biDir[n] *= abs(((1-s_b) + s_b * np.dot(biDir[n], micDir)))
    varDir[n] *= abs(((1-s_v) + s_v * np.dot(varDir[n], micDir)))

fig = plt.figure()
ax = fig.gca()
ax.axhline(y=0, color='k', linewidth=.5)
ax.axvline(x=0, color='k', linewidth=.5)
ax.set_xticks(np.arange(-1, 1.1, 0.2))
ax.set_yticks(np.arange(-1, 1.1, 0.2))


def render(signal):
    plt.plot(*signal.T)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.minorticks_on()
    plt.grid(True, linestyle='--')

    #plt.show()


# plots:
render(omniDir)
render(cardioid)
render(biDir)
render(varDir)

plt.show()