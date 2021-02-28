import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.fftpack import fft, fftfreq
from IPython.lib.display import Audio
plt.style.use("ggplot")


file = "E:/OneDrive/Desktop/wall zeros.wav"
[f_s, track] = wavfile.read(file)
# just left channel:
data = track[:, 0]


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
            max = sig[i]
    return sig/max


def createSinusSignal(frq):
    # Number of sample points
    fs = 44100
    N = fs
    # sample spacing
    T = 1.0 / fs
    myX = np.linspace(0.0, N * T, N, endpoint=False)
    myY = 0.3*np.sin(frq * 2.0 * np.pi * myX)
    return myY


# data = createSinusSignal(300)
# f_s = 44100

# plot im Zeitbereich:
t = np.linspace(0, data.size/f_s, data.size)
x_t = normalize(data)

fig = plt.figure()
ax = fig.gca()
# ax.axvline(x=0.1, color='b', linewidth=1.5, ls='--')
# plt.text(0.1 - 0.04, 0.11, 'Beamtracing', color='b', rotation='0', size='16')
# plt.text(0.11, 0.11, 'Raytracing', color='b', rotation='0', size='16')

plt.plot(t, data)
plt.minorticks_on()
# plt.title("Data im Zeitbereich")
plt.xlim(0.02675, 0.0275)
plt.xlabel("t in s")
plt.ylabel("$h(t)$")
plt.show()

# spectrum
o = 8192*4
ohalf = int(o/2)
N = len(data)
y = abs(np.fft.fft(data, o))
yl = 20*np.log10(normalize(y))
print(len(yl))
x = makeFAxis(yl, f_s)
plt.plot(x[:ohalf], yl[:ohalf], linewidth=0.5)
plt.title("Frequenzbereich bis $f_s/ 2$")
plt.xlabel("$frq$")
plt.ylim(-30, 5)
plt.xlim([10, f_s/2])
plt.semilogx()

plt.fill_between(x, -70, yl, alpha=0.5)
plt.fill_between(x, 5, yl, alpha=0.7, facecolor="black")

plt.show()
