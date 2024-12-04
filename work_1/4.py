import numpy as np
import matplotlib.pyplot as plt

def qpsk_modulate(data_bits):
    symbols = []
    for i in range(0, len(data_bits), 2):
        bit_pair = data_bits[i : i + 2]
        if np.array_equal(bit_pair, [0, 0]):    symbols.append(np.exp(1j * 0))       # 0 radians
        
        elif np.array_equal(bit_pair, [0, 1]):  symbols.append(np.exp(1j * (np.pi / 2)))  # pi/2 radians
        
        elif np.array_equal(bit_pair, [1, 0]):  symbols.append(np.exp(1j * np.pi))   # pi radians
        
        elif np.array_equal(bit_pair, [1, 1]):  symbols.append(np.exp(1j * (3 * np.pi / 2)))  # 3pi/2 radians
        
    return np.array(symbols)

def AWGN(signal, signal_to_noise_ratio):
    target_noise_watts = 10 ** (signal_to_noise_ratio / 10)
    signal_power = np.mean(np.abs(signal) ** 2)
    noise_power = signal_power / target_noise_watts
    noise = np.sqrt(noise_power) * (np.random.normal(size = signal.shape) + 1j * np.random.normal(size = signal.shape))
    return signal + noise

N = 8 # amount of subcarriers
symbols_count = 10
bits_count = N * symbols_count * 2 #QPSK => *2
signal_to_noise_ratio = 25
data_bits = np.random.randint(0, 2, bits_count)

# Модуляция QPSK
qpsk_symbols = qpsk_modulate(data_bits)  * np.exp(1j * np.pi / 4) 
ifft_symbols = np.fft.ifft(qpsk_symbols)
awgn_symbols = AWGN(ifft_symbols, signal_to_noise_ratio)
received = np.fft.fft(awgn_symbols)


plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(np.real(ifft_symbols), label = 'Re')
plt.plot(np.imag(ifft_symbols), label = 'Im')
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.scatter(np.real(received), np.imag(received))
plt.grid()
plt.show()     
