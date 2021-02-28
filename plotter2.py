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
            max = sig[i]
    return sig/max

# https://www.schweizer-fn.de/stoff/akustik/absorptionsfaktoren.php


f1 = 125
f2 = 250
f3 = 500
f4 = 1000
f5 = 2000
f6 = 4000

a1 = 0.9  # in Band 0-200
a2 = 0.8  # in Band 200-400
a3 = 1.0  # in Band 400-800
a4 = 0.3  # in Band 800-1600
a5 = 0.3  # in Band 1600-3200
a6 = 0.5  # in Band 3200 - fs/2

fs = 24000
o = 2048*2
n = np.linspace(0, 1, fs)
# y1 = a1*np.sin(2*np.pi*f1*n)
y = a1*np.sin(2*np.pi*f1 * n) + a2*np.sin(2*np.pi*f2 * n) + a3*np.sin(2*np.pi*f3 * n) + a4*np.sin(2*np.pi*f4 * n) + a5*np.sin(2*np.pi*f5 * n) + a6*np.sin(2*np.pi*f6 * n)

yint = 0

# Interpolation:
Yint = np.zeros(2*o)
for k in range(0, o):
    frq = k * fs/o
    if frq < f1:
        deltacoeff = a2 - a1
        bandw = f2 - f1
        deltafband = frq - f1
        a = a1 + ((deltafband * deltacoeff)/bandw)
    elif f1 <= frq < f2:
        deltacoeff = a2 - a1
        bandw = f2 - f1
        deltafband = frq - f2
        a = a2 + ((deltafband * deltacoeff) / bandw)
    elif f2 <= frq < f3:
        deltacoeff = a3 - a2
        bandw = f3 - f2
        deltafband = frq - f3
        a = a3 + ((deltafband * deltacoeff) / bandw)
    elif f3 <= frq < f4:
        deltacoeff = a4 - a3
        bandw = f4 - f3
        deltafband = frq - f4
        a = a4 + ((deltafband * deltacoeff) / bandw)
    elif f4 <= frq < f5:
        deltacoeff = a5 - a4
        bandw = f5 - f4
        deltafband = frq - f5
        a = a5 + ((deltafband * deltacoeff) / bandw)
    elif f5 <= frq < f6:
        deltacoeff = a6 - a5
        bandw = f6 - f5
        deltafband = frq - f6
        a = a6 + ((deltafband * deltacoeff) / bandw)
    else:
        deltacoeff = 0 - a6
        bandw = fs - f6
        deltafband = frq - fs
        a = 0 + ((deltafband * deltacoeff) / bandw)
    # zeitbereich:
    yy = a*np.sin(2*np.pi*frq*n)
    yint += yy
    # Frq-Bereich:
    Yint[k] = a
    if k > 0:
        # Spiegelspektrum
        Yint[len(Yint)-k] = a


# plot

# spektrum der 6 Ausgangsfrequenzen:
Y = normalize(abs(np.fft.fft(y, o)))
YAxis = makeFAxis(Y, fs)
plt.stem(YAxis[0:int(o/2)], Y[0:int(o/2)])
plt.semilogx()
plt.title("Spektrum der 6 Ausgangsfrequenzen")
plt.show()

# Zeitbereich:
T1 = int(fs/f1)
# plt.plot(n[0:T1], normalize(y1[0:T1]))
plt.plot(n[0:T1], normalize(y[0:T1]))
plt.title("Addition der 6 Frequenzen im Zeitbereich")
plt.show()

T1 = int(fs/f1)
# plt.plot(n[0:T1], normalize(y1[0:T1]))
plt.plot(n[0:T1], normalize(yint[0:T1]))
plt.title("Addition aller interpolierten Frequenzen")
plt.show()

# Spektrum von yint
Y2 = normalize(abs(np.fft.fft(yint, o)))
Y2Axis = makeFAxis(Y2, fs)
plt.stem(Y2Axis[0:int(o/2)], Y2[0:int(o/2)])
plt.semilogx()
plt.title("Interpolation im Zeitbereich und FFT")
plt.show()

# interpolation im Frequenzbereich:
plt.stem(np.linspace(0, fs, 2*o), Yint)
plt.semilogx()
plt.title("Interpolation im Frequenzbereich (Wie in BA)")
plt.xlim(10, fs/2)
plt.show()

# ifft von dem selbst interpolierten Spektrum:
ir = np.fft.ifft(Yint.copy())
nn = np.linspace(0, (2*o)/fs, 2*o)
plt.plot(nn, ir.real)
plt.title("IFFT des interpolierten Spektrums ohne shift:")
plt.show()

# shift:
shifted = np.zeros(2*o, complex)
shifted[:o] = ir[o:]
shifted[o:] = ir[:o]
plt.plot(nn, shifted.real)
plt.title("mit shift Zeitbereich")
lower = ((len(shifted)/2)-128)/fs
upper = ((len(shifted)/2)+128)/fs
plt.xlim(lower, upper)
plt.show()

# fft von der ifft als Vergleich (ohne shift):
irfft = abs(np.fft.fft(ir, 2*o))
nfft = np.linspace(0, fs, 2*o)
plt.stem(nfft[0:int(o)], irfft[0:int(o)])
plt.title("Spektrum impuls ohne shift")
plt.semilogx()
plt.show()

# fft von der ifft als Vergleich (mit shift):
irfftshifted = abs(np.fft.fft(shifted, 2*o))
plt.stem(nfft[0:int(o)], irfftshifted[0:int(o)])
plt.title("Spektrum impuls mit shift")
plt.semilogx()
plt.show()

# Butterworth TP/HP:
spectrum = Yint.copy()
order = 4
fGrenzO = 3000
fGrenzU = 60
HBWTP = np.array([])
HBWHP = np.array([])
for fIndex in range(1, 2*o):
    fval = fIndex * fs/(2*o)
    tpbw = 1 / (1 + pow(fval/fGrenzO, 2*order))
    hpbw = 1 / (1 + pow(fGrenzU/fval, 2*order))
    spectrum[fIndex] *= tpbw
    spectrum[fIndex] *= hpbw
    HBWTP = np.append(HBWTP, tpbw)
    HBWHP = np.append(HBWHP, hpbw)

#plt.stem(nfft[0:int(o)], spectrum[0:int(o)])
plt.plot(nfft[0:int(o)], HBWTP[0:int(o)])
plt.plot(nfft[0:int(o)], HBWHP[0:int(o)])
plt.title("Butterworth HP & TP mit n = 4")
plt.semilogx()
plt.show()

hbwtp = abs(np.fft.fft(HBWTP, o))
#hbwhp = np.fft.fft(HBWTP, o)
nnfft = np.linspace(0, fs, o)
plt.plot(nnfft, hbwtp)
#plt.plot(nnfft, hbwhp)
plt.xlim(0, 256)
plt.show()
