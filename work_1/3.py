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
delay = 1
t = np.linspace(-1*T, 4*T, 512)

sequence_video = np.zeros_like(t)
for i in range(4):
    start = i * (T + delay)
    sequence_video += np.where((t >= start) & (t < start + delay), A, 0)

freqs_seq_video = np.arange(len(sequence_video))

spectr_seq_video = DFT(sequence_video)

plt.figure(3, label = 'Sequence of videoimpulses')
plt.subplot(2, 1, 1)
plt.plot(t, sequence_video)
plt.grid(True)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal")
plt.show()


plt.subplot(2, 1, 2)
plt.plot( freqs_seq_video, np.abs(spectr_seq_video))
plt.grid(True)
plt.xlabel("Frequencies")
plt.ylabel("Amplitude")
plt.title("Spectr")
plt.show()