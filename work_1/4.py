import numpy as np
import matplotlib.pyplot as plt

N = 64 # amount of subcarriers
symbols_count = 5
bits_count = N * symbols_count * 2 #QPSK => *2
signal_to_noise_ratio = 15
bits = np.random.randint(0, 2, bits_count)

symbols = np.zeros((symbols_count, N), dtype = complex)

for i in range(symbols_count):
    for j in range(0, N, 2): #step = 2 since it's QPSK
        bit1 = bits[i * 2 * N + j]
        bit2 = bits[i * 2 * N + j + 1]
        
        if (bit1, bit2) == (0, 0): symbols[i, j // 2] = 1 + 1j #00

        if (bit1, bit2) == (0, 1): symbols[i, j // 2] = -1 + 1j #01
        
        if (bit1, bit2) == (1, 1): symbols[i, j // 2] = -1 - 1j #11
        
        if (bit1, bit2) == (1, 0): symbols[i, j // 2] = 1 - 1j #10

ifft_symbols = np.fft.ifft(symbols)

def AWGN(signal, signal_to_noise_ratio):
    target_noise_watts = 10 ** (signal_to_noise_ratio / 10)
    signal_power = np.mean(np.abs(signal) ** 2)
    noise_power = signal_power / target_noise_watts
    noise = np.sqrt(noise_power) * (np.random.normal(size = signal.shape) + 1j * np.random.normal(size = signal.shape))
    return signal + noise

ifft_awgn_symbols = AWGN(ifft_symbols, signal_to_noise_ratio)

received = np.fft.fft(ifft_awgn_symbols)

plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(np.real(ifft_symbols[0]), label = 'Re')
plt.plot(np.imag(ifft_symbols[0]), label = 'Im')
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.scatter(np.real(received[0]), np.imag(received[0]))
plt.grid()
plt.show()       
    