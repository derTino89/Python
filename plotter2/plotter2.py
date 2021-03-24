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


f1 = 125
f2 = 250
f3 = 500
f4 = 1000
f5 = 2000
f6 = 4000

a1 = 0.9
a2 = 0.8
a3 = 1.0
a4 = 0.3
a5 = 0.3
a6 = 0.5

fs = 48000
o = int(fs/f1)
ohalf = int(o/2)
order = 4
fGrenzO = 12000
fGrenzU = 100
HBWBP = np.array([])
n = np.linspace(0, 1, fs)
y = a1*np.sin(2*np.pi*f1 * n) + a2*np.sin(2*np.pi*f2 * n) + a3*np.sin(2*np.pi*f3 * n) + a4*np.sin(2*np.pi*f4 * n) + a5*np.sin(2*np.pi*f5 * n) + a6*np.sin(2*np.pi*f6 * n)
nfft = np.linspace(0, fs, o)
yint = 0

# Interpolation:
Hraw = np.empty(o, dtype=complex)
HrawNoShift = np.empty(o, dtype=complex)
Hfilter = np.empty(o, dtype=complex)

for k in range(0, ohalf+1):
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
        deltacoeff = 0
        bandw = fs/2 - f6
        deltafband = frq - fs/2
        a = a6 + ((deltafband * deltacoeff) / bandw)
    # zeitbereich:
    yy = a*np.sin(2*np.pi*frq*n)
    yint += yy
    # Frq-Bereich:
    # Filter:
    if k == 0:
        tpbw = 0
        hpbw = 0
    else:
        tpbw = 1 / np.sqrt((1 + pow(frq / fGrenzO, 2 * order)))
        hpbw = 1 / np.sqrt((1 + pow(fGrenzU / frq, 2 * order)))
    # gain:
    # shift impuls to mid of array by adding phase with multiplying the abs() with exp(2pi k*m/N), where k is the frequenzstützstelle n, N is amount of stützstellen or blocksize, m is the amount we want to shift
    HBWBP = np.append(HBWBP, (tpbw * hpbw))
    Hraw[k] = a * np.exp(-2j * np.pi * (k * ohalf) / o)
    HrawNoShift[k] = a
    Hfilter[k] = a * np.exp(-2j * np.pi * (k * ohalf) / o) * tpbw * hpbw
    if k < o:
        # Spiegelspektrum
        Hraw[-k] = a * np.exp(2j * np.pi * (k * ohalf) / o)
        HrawNoShift[-k] = a
        Hfilter[-k] = a * np.exp(2j * np.pi * (k * ohalf) / o) * tpbw * hpbw
# plot
# Zeitbereich:
T1 = int(fs/f1)
plt.subplot(1, 2, 1)
plt.plot(n[0:T1], normalize(y[0:T1]))
#plt.title("Addition der 6 Frequenzen im Zeitbereich")
plt.xlabel("$t$")
plt.ylabel("$y(t)$")
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')

# spektrum der 6 Ausgangsfrequenzen:
Y = normalize(abs(np.fft.fft(y, o)))
YAxis = makeFAxis(Y, fs)
plt.subplot(1, 2, 2)
plt.stem(YAxis[0:ohalf], Y[0:ohalf])
plt.xlabel('$\omega$')
plt.ylabel("$Y(j\omega$)")
plt.semilogx()
plt.xlim(10, fs/2)
plt.grid(True, linestyle='--', linewidth='0.3')
#plt.title("Spektrum der 6 Ausgangsfrequenzen")
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')
plt.show()


# interpolation im Frequenzbereich:
fig = plt.figure()
ax = fig.gca()
specAxis = np.linspace(0, fs, o)
plt.stem(specAxis, abs(Hraw))
plt.semilogx()
plt.xlim(10, fs/2)
plt.grid(True, linestyle='--', linewidth='0.3')
plt.xlabel('$\omega$')
plt.ylabel("$H(j\omega$)")
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')
plt.show()


# Butterworth TP/HP:
fig2 = plt.figure()
gs = fig2.add_gridspec(2, 2)
f2_1 = fig2.add_subplot(gs[1, :])
plt.stem(abs(Hfilter[0:ohalf]))
plt.plot(HBWBP[0:ohalf], color='r')
plt.grid(True, linestyle='--', linewidth='0.3')
plt.xlabel('$k$')
plt.ylabel("$H[k]$")
#plt.semilogx()
#plt.xlim(10, fs/2)
f2_2 = fig2.add_subplot(gs[0, 0])
#ax2 = fig2.gca()
#ax2.axvline(x=fGrenzU, color='r', linestyle='--')
#ax2.axvline(x=fGrenzO, color='r', linestyle='--')
plt.plot(nfft[0:ohalf], HBWBP[0:ohalf], color='r')
plt.grid(True, linestyle='--', linewidth='0.3')
plt.xlabel('$\omega$')
plt.ylabel("$H_{BWBP}(j\omega$)")
#plt.semilogx()
#plt.xlim(10, fs/2)
f2_3 = fig2.add_subplot(gs[0, 1])
plt.stem(nfft[0:ohalf], abs(Hraw[0:ohalf]))
plt.grid(True, linestyle='--', linewidth='0.3')
plt.xlabel('$\omega$')
plt.ylabel("$H(j\omega$)")
#plt.semilogx()
#plt.xlim(10, fs/2)
plt.show()

plt.subplot(2, 1, 1)
plt.stem(np.abs(Hfilter))
plt.xlabel('$k$')
plt.ylabel("$H[k]$")
#ax[1].stem(np.unwrap(np.angle(Hfilter)))
plt.subplot(2, 1, 2)
plt.stem(np.angle(HrawNoShift))
plt.xlabel('$k$')
plt.ylabel("$\phi[k]$")
plt.ylim(-1., 1.)
plt.show()

plt.subplot(2, 1, 1)
plt.stem(np.abs(Hfilter))
plt.xlabel('$k$')
plt.ylabel("$H[k]$")
#ax[1].stem(np.unwrap(np.angle(Hfilter)))
plt.subplot(2, 1, 2)
plt.stem(np.angle(Hfilter))
plt.xlabel('$k$')
plt.ylabel("$\phi[k]$")
plt.show()


# ifft von dem selbst interpolierten Spektrum ohne shift:
ir = np.fft.ifft(HrawNoShift.copy())
plt.subplot(1, 2, 1)
plt.plot(np.real(ir))
plt.xlabel("n")
plt.ylabel("h[n]")
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')


# ifft von dem selbst interpolierten Spektrum mit shift:
ir = np.fft.ifft(Hraw.copy())
plt.subplot(1, 2, 2)
plt.plot(np.real(ir))
plt.xlabel("n")
plt.ylabel("h[n - N/2]")
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')
plt.show()

# window:
window = np.empty(o)
for nIndex in range(0, o):
    window[nIndex] = 0.35875 - (0.48829 * (np.cos(2 * np.pi * nIndex / o))) + (0.14128 * (np.cos(4 * np.pi * nIndex / o))) - (0.01168 * (np.cos(6 * np.pi * nIndex / o)))
plt.subplot(1, 2, 1)
plt.plot(window, 'r')
plt.xlabel("n")
plt.ylabel("hw[n]")

ir = np.fft.ifft(Hfilter.copy())
windowedIR = ir * window
plt.subplot(1, 2, 2)
plt.plot(ir)
plt.plot(windowedIR)
plt.xlabel("n")
plt.ylabel("h[n-N/2]")
plt.show()


# fft von der ifft als Vergleich:
irfft = abs(np.fft.fft(ir, o))
plt.stem(nfft, irfft)
plt.title("Spektrum impuls")
plt.semilogx()
plt.xlim(10, fs/2)
plt.grid(True, linestyle='--', linewidth='0.3')
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')
plt.show()
# hbwbp = abs(np.fft.fft(HBWBP, o))
# nnfft = np.linspace(0, fs, o)
# plt.plot(nnfft, hbwbp)
# plt.xlim(0, 1024)
# mng = plt.get_current_fig_manager()
# mng.window.state('zoomed')
# plt.show()
