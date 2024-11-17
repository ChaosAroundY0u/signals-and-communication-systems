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
t = np.linspace(-1*T, 4*T, 512)

signal_video = np.where((t >= 0) & (t <= 0 + T), A, 0)

freqs_video = np.arange(len(signal_video))

spectr_video = DFT(signal_video)

plt.figure(1, label = "Videoimpuls")
plt.subplot(2, 1, 1)
plt.plot(t, signal_video)
plt.grid(True)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal")
plt.show()

plt.subplot(2, 1, 2)
plt.plot(freqs_video, np.abs(spectr_video))
plt.grid(True)
plt.xlabel("Frequencies")
plt.ylabel("Amplitude")
plt.title("Spectr")
plt.show()