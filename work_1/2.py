import numpy as np
import matplotlib.pyplot as plt

def DFT(signal):
    N = len(signal)
    ans = np.zeros(N, dtype = complex)
    for k in range (N):
        for n in range (N):
            ans[k] += signal[n] * np.exp(-1j * 2*np.pi * k * n / N)
    return ans

T = 5
A = 2
f0 = 2
t = np.linspace(-1*T, 4*T, 512)


signal_radio = np.where((t >= 0) & (t <= 0 + T), np.cos(f0*t), 0)

freqs_radio = np.arange(len(signal_radio))

spectr_radio = DFT(signal_radio)

plt.figure(2, label = 'Radioimpuls')
plt.subplot(2, 1, 1)
plt.plot(t, signal_radio)
plt.grid(True)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal")
plt.show()

plt.subplot(2, 1, 2)
plt.plot(freqs_radio, np.abs(spectr_radio))
plt.grid(True)
plt.xlabel("Frequencies")
plt.ylabel("Amplitude")
plt.title("Spectr")
plt.show()