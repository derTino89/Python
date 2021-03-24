import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
plt.style.use("ggplot")


fileA = "E:/OneDrive/Desktop/BA RIR/ba_wall_RT_6mio.wav"
[f_sA, trackA] = wavfile.read(fileA)
# just left channel:
dataA = trackA[:, 0]

#fileB = "E:/OneDrive/Desktop/BA RIR/ba_wall_RT_125k.wav"
fileB = "E:/OneDrive/Desktop/BA RIR/class test 125k 4 no direct sound.wav"
[f_sB, trackB] = wavfile.read(fileB)
# just left channel:
dataB = trackB[:, 0]
dataBR = trackB[:, 1]

fileC = "E:/OneDrive/Desktop/BA RIR/ba_wall_BT.wav"
[f_sC, trackC] = wavfile.read(fileC)
# just left channel:
dataC = trackC[:, 1]

# rechter kanal:
#dataB = trackA[:, 1]

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
tA = np.linspace(0, (dataA.size-1)/f_sA, dataA.size)
tB = np.linspace(0, (dataB.size-1)/f_sB, dataB.size)
tC = np.linspace(0, (dataC.size-1)/f_sC, dataC.size)


# # singleplot (fileB):
# plt.plot(tB, dataB)
# plt.show()

# stereoplot (fileB):
plt.subplot(2,1,1)
plt.plot(tB, dataB)
plt.minorticks_on()
plt.xlim(0, 0.7)
plt.xlabel("$t$")
plt.ylabel("$h_L(t)$")

plt.subplot(2,1,2)
plt.plot(tB, dataBR)
plt.minorticks_on()
plt.xlim(0, 0.7)
plt.xlabel("t")
plt.ylabel("$h_R(t)$")

plt.show()

# # Multiplot:
# # file A:
# plt.subplot(1,3,1)
# plt.stem(tA, dataA)
# plt.plot(tA, dataA, 'b', linewidth=0.5)
# plt.minorticks_on()
# # wall-Szene:
# plt.xlim(0.1942, 0.1948)
# # direct sound:
# #plt.xlim(0.029, 0.0293)
# # class-Szene:
# #plt.xlim(0, 1)
# plt.ylim(-0.001, 0.009)
# plt.xlabel("t in s")
# plt.ylabel("$h(t)$")
# plt.title("RT: $N = 6.000.000$")
#
# # file B:
# plt.subplot(1,3,2)
# plt.stem(tB, dataB)
# plt.plot(tB, dataB, 'b', linewidth=0.5)
# #plt.stem(dataB)
# # plt.plot(dataB, 'b', linewidth=0.5)
# plt.minorticks_on()
# # wall-Szene:
# plt.xlim(0.1942, 0.1948)
# # direct sound:
# #plt.xlim(0.029, 0.0293)
# # class-Szene:
# #plt.xlim(0, 1)
# plt.xlabel("t in s")
# plt.ylabel("$h(t)$")
# plt.ylim(-0.001, 0.009)
# plt.title("RT: $N = 125.000$")
#
# # file C:
# plt.subplot(1,3,3)
# plt.stem(tC, dataC)
# plt.plot(tC, dataC, 'b', linewidth=0.5)
# #plt.stem(dataC)
# # plt.plot(dataC, 'b', linewidth=0.5)
# plt.minorticks_on()
# # wall-Szene:
# plt.xlim(0.1942, 0.1948)
# # direct sound:
# #plt.xlim(0.029, 0.0293)
# # class-Szene:
# #plt.xlim(0, 1)
# plt.xlabel("t in s")
# plt.ylabel("$h(t)$")
# plt.ylim(-0.001, 0.009)
# plt.title("BT: 1. Ordnung")
#
# plt.show()



# spectrum
o = 2048*4
ohalf = int(o/2)
N = len(dataB)
y = abs(np.fft.fft(dataB, o))
yl = 20*np.log10(y)
#print(len(yl))
x = makeFAxis(yl, f_sB)
plt.plot(x[:ohalf], yl[:ohalf], linewidth=0.5)
#plt.title("Frequenzbereich bis $f_s/ 2$")
plt.xlabel("$\omega / 2\pi$")
plt.ylabel("$H(j\omega)$")
plt.ylim(-60, 5)
plt.xlim([10, f_sB/2])
plt.semilogx()

plt.fill_between(x, -70, yl, alpha=0.5)
#plt.fill_between(x, 5, yl, alpha=0.7, facecolor="black")

plt.show()
