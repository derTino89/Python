import numpy as np
import matplotlib.pyplot as plt


def makeFAxis(spectrum, Fsample):
    mySize = spectrum.size
    fk_Axis = np.zeros(mySize)
    for k in range(mySize):
        fk_Axis[k] = (k/mySize)*Fsample
    return fk_Axis


def normalize(sig):
    max = 0
    for i in range(0, sig.size):
        if sig[i] > max:
            max = \
                sig[i]
    return sig/max


def complex_stem(H):
    fig, ax = plt.subplots(1, 2)
    ax[0].stem(np.abs(H))
    ax[1].stem(np.angle(H))
    plt.show()

# https://www.schweizer-fn.de/stoff/akustik/absorptionsfaktoren.php




fs = 8000

fGrenzO = 500
fGrenzU = 100
N = 1024

t = np.linspace(-2, 2, N)
d = np.zeros(N)
d[int(N/2)] = 1
D = np.fft.fft(d, N)

w = np.ones(N)
TP = np.ones(N)
for k in range(0,N):
    w[k] = k*fs/N
    TP[k] *= 1 / np.sqrt((1 + pow(w[k] / fGrenzO, 2 * N))) * np.exp(-1j * 2 * np.pi * N/2 * k / N)
plt.subplot(2,2,1)
plt.plot(t, d)
plt.xlabel('t')
plt.ylabel('$x(t)$')
plt.subplot(2,2,2)
plt.plot(w,abs(D))
plt.semilogx()
plt.ylim(0,)
plt.xlabel('$\omega$')
plt.ylabel('$X(j\omega)$')
plt.subplot(2,2,4)
plt.semilogx()
plt.plot(w,abs(TP))
tp = np.fft.ifft(TP)
plt.xlabel('$\omega$')
plt.ylabel('$H(j\omega)$')
plt.subplot(2,2,3)
plt.plot(t, np.real(tp))
plt.xlabel('$t$')
plt.ylabel('$h(t)$')
plt.show()


