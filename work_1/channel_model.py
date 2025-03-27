import numpy as np
import matplotlib.pyplot as plt

def qpsk_modulate(data_bits):
    bit1 = data_bits[0::2]
    bit2 = data_bits[1::2]
    x = 1/np.sqrt(2) * (1 - 2*bit1)
    y = 1/np.sqrt(2) * (1 - 2*bit2)
    qpsk_symbols = x + 1j * y

    return qpsk_symbols

def AWGN(signal, signal_to_noise_ratio):
    noise_power = 10 ** (-signal_to_noise_ratio / 10)
    print(noise_power)
    noise = np.sqrt(noise_power / 2) * (np.random.normal(size=signal.shape) + 1j * np.random.normal(size=signal.shape))
    return signal + noise

def tworay_channel(signal, delay_samples, alpha):
    
    fs = 32.76e6
    delay_smpls = int(30 * fs / 10e9)
    signal = signal / np.sqrt(1 + alpha**2)
    delayed_signal = np.concatenate((np.zeros(delay_samples), signal[:-delay_samples])) * alpha / np.sqrt(1 + alpha**2)
    summ = signal + delayed_signal
    return summ


# Параметры
N = 512  # количество поднесущих
symbols_count = 10
cp_len = 16
bits_count = N * symbols_count * 2  # QPSK => *2
signal_to_noise_ratio = 15



data_bits = np.random.randint(0, 2, bits_count)

# QPSK
qpsk_symbols = qpsk_modulate(data_bits)

pilot_symbol = qpsk_symbols[:N] # первый символ из 10 - пилотный

qpsk_symbols = qpsk_symbols.reshape(symbols_count, N)

# IFFT
ifft_symbols = np.fft.ifft(qpsk_symbols, norm="ortho")

# + циклический префикс
cp_ifft_symbols = np.hstack((ifft_symbols[...,-cp_len:], ifft_symbols))


delay_sampl = 10
tworay_signal = tworay_channel(cp_ifft_symbols.ravel(), delay_sampl, 0.5)


noise_power = 10 ** (-signal_to_noise_ratio / 10)

# + шум
awgn_symbols = AWGN(tworay_signal, signal_to_noise_ratio)

# - циклический префикс
no_cp_symbols = awgn_symbols.reshape(symbols_count, -1)[..., cp_len:]


received = np.fft.fft(no_cp_symbols, norm="ortho")


h_array = (received[0] / pilot_symbol)
W_mmse = np.conj(h_array) / (h_array * np.conj(h_array) + (noise_power))
print(W_mmse)
#h_array = np.dot(np.dot(np.linalg.inv(np.dot(pilot_symbol.T, pilot_symbol)), pilot_symbol.T), received[0])
equaliz_received = received * W_mmse
equaliz2 = received / h_array

plt.figure(2)
plt.scatter(np.real(received), np.imag(received))

plt.scatter(np.real(equaliz_received), np.imag(equaliz_received))
plt.scatter(np.real(equaliz2), np.imag(equaliz2))
plt.grid()
plt.show()
plt.figure(3)
plt.scatter(np.real(equaliz_received), np.imag(equaliz_received))
plt.grid()
plt.show()
