import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

def qpsk_modulate(data_bits):
    return 1/np.sqrt(2) * ((1-2*data_bits[0::2]) + 1j*(1-2*data_bits[1::2])) #TS 38.211 - 5.1.3

def AWGN(signal, signal_to_noise_ratio):
    noise_power = 10 ** (-signal_to_noise_ratio / 10)
    print(noise_power)
    noise = np.sqrt(noise_power / 2) * (np.random.normal(size=signal.shape) + 1j * np.random.normal(size=signal.shape))
    return signal + noise

def Normal_Distribution(z):
    return (1 / np.sqrt(2*np.pi)) * np.exp(-0.5 * z ** 2)

def TDLA(signal):
    #tdl-a parameters
    delays_ns = np.array([0, 38.19, 40.25, 58.68, 46.10, 53.75, 67.08]) #nominal delay spread

    powers_db = np.array([-13, 0, -2.2, -4, -6, -8.2, -9.9])
    #from db to linear
    powers_linear = 10 ** (powers_db / 10)

    tdl_powers = powers_linear
    tdl_powers /= np.max(tdl_powers)
    #impulse response
    sampling_rate = 2048 * 15e3 #N_fft * delta_f = 30 720 000
    delays_samples = np.array([np.round(delays_ns[i] * sampling_rate / 1e9).astype(int) for i in range(len(delays_ns))])
    
    impulse_response = np.zeros_like(powers_linear)
    a = np.array(list(zip(delays_samples, powers_linear)))
    print(a)
    #for delay, power in zip(delays_samples, powers_linear):
    #   if delay < len(impulse_response):
    for i in range(len(powers_linear)):
            impulse_response[i] +=  np.sqrt(powers_linear[i]) * Normal_Distribution(powers_linear[i])
    #norm
    output_signal = np.zeros_like(signal)  
    norm_koef = np.sqrt(np.sum(np.abs(impulse_response))**2)
    normed_ir = impulse_response / norm_koef
    for i, koef in enumerate(normed_ir):
        if i == 0:
            output_signal += signal * koef
        else:
            delayed_signal = np.concatenate((np.zeros(i), signal[:-i])) * koef
            output_signal += delayed_signal
    return output_signal

#ввести частоту семплирования
sampling_rate_Hz = 15e3 #delta f_ref
"""mu = 0 - subcarrier spacing configutation"""
#расстояние между поднесущими 15кгц

# Параметры
N = 2048  # количество поднесущих
symbols_count = 10
cp_len_first = 160
cp_len = 144
bits_count = N * symbols_count * 2  # QPSK => *2
signal_to_noise_ratio = 100

data_bits = np.random.randint(0, 2, bits_count)

# QPSK
qpsk_symbols = qpsk_modulate(data_bits)

pilot_symbol = qpsk_symbols[:N] # первый символ из 10 - пилотный

qpsk_symbols = qpsk_symbols.reshape(symbols_count, N)

# IFFT
ifft_symbols = np.fft.ifft(qpsk_symbols, norm="ortho")
print(np.hstack((ifft_symbols[0, -cp_len_first:], ifft_symbols[0])))
cp_ifft_symbols = np.array([])
# + циклический префикс
for i in range(symbols_count):
    if i == 0:
        cp_ifft_symbols = np.append(np.hstack((ifft_symbols[i, -cp_len_first:], ifft_symbols[i])), cp_ifft_symbols)
    else:
        cp_ifft_symbols = np.append(np.hstack((ifft_symbols[i,-cp_len:], ifft_symbols[i])), cp_ifft_symbols)


TDL_A_signal = TDLA(cp_ifft_symbols)


noise_power = 10 ** (-signal_to_noise_ratio / 10)

# + шум
awgn_symbols = AWGN(TDL_A_signal, signal_to_noise_ratio)

#- циклический префикс
no_cp_symbols_first = (awgn_symbols[:(N+cp_len_first)])[cp_len_first:]
awgn_symbols = awgn_symbols[(N+cp_len_first):]
no_cp_symbols = awgn_symbols.reshape(symbols_count-1, -1)[..., cp_len:]
no_cp_symbols = np.vstack((no_cp_symbols_first, no_cp_symbols))

#no_cp_symbols = awgn_symbols.reshape(symbols_count, -1)[..., cp_len:]


received = np.fft.fft(no_cp_symbols, norm="ortho")


h_array = (received[0] / pilot_symbol)
W_mmse = np.conj(h_array) / (h_array * np.conj(h_array) + (noise_power))

equaliz_received = received * W_mmse

plt.figure(2)
plt.scatter(np.real(received), np.imag(received))

plt.scatter(np.real(equaliz_received), np.imag(equaliz_received))

plt.grid()
plt.show()
