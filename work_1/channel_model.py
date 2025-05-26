import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

def qpsk_modulate(data_bits):
    return 1/np.sqrt(2) * ((1-2*data_bits[0::2]) + 1j*(1-2*data_bits[1::2])) #TS 38.211 - 5.1.3

def AWGN(signal, noise_power):
    print("noise power:",noise_power)
    noise = np.sqrt(noise_power / 2) * (np.random.normal(size=signal.shape) + 1j * np.random.normal(size=signal.shape))
    return signal + noise
    
def TDLa(signal, delay_spread, N_fft, delta_f):
    
    delays_ns = np.array([0, .3819, .4025, .5868, .4610, .5375, .6708]) * delay_spread #nominal delay spread
    powers_db = np.array([-13, 0, -2.2, -4, -6, -8.2, -9.9])
    powers_linear = 10 ** (powers_db / 10)

    tdl_powers = powers_linear
    tdl_powers /= np.max(tdl_powers)

    sampling_rate = N_fft * delta_f #N_fft * delta_f = 30 720 000
    delays_samples = np.array([np.round(delays_ns[i] * sampling_rate / 1e9).astype(int) for i in range(len(delays_ns))])
    
    impulse_response = np.zeros(np.max(delays_samples) + 1, dtype=complex)

    for delay, power in zip(delays_samples, powers_linear):
        impulse_response[delay] += np.sqrt(power) * (np.random.normal(0, 1/np.sqrt(2)) + 1j * np.random.normal(0, 1/np.sqrt(2)))   

    output_signal = np.convolve(signal, impulse_response, mode = "same")

    print("impulse response:", impulse_response)   
    return output_signal, powers_linear, delays_samples

def covariance_matrix(N_fft, T_cp, T_OFDM):
    alpha = T_cp / T_OFDM
    R_hh = np.zeros((N_fft, N_fft), dtype=complex)
    for k in range(N_fft):
        for m in range(N_fft):
            if k==m: R_hh[k, m] = 1.0
        else:
            delta = k-m
            phase = -np.pi * delta * alpha
            R_hh[k, m] = (1 / alpha) * np.sinc(delta * alpha) * np.exp(1j * phase)
    return R_hh

# Параметры
N = 512  # количество поднесущих
symbols_count = 10
cp_len = 64
bits_count = N * symbols_count * 2  # QPSK => *2
signal_to_noise_ratio = 15
delay_spread = 1
delta_f = 15e3 # Hz

data_bits = np.random.randint(0, 2, bits_count)

# QPSK
qpsk_symbols = qpsk_modulate(data_bits)

pilot_symbol = qpsk_symbols[:N] # первый символ из 10 - пилотный

qpsk_symbols = qpsk_symbols.reshape(symbols_count, N)

# IFFT
ifft_symbols = np.fft.ifft(qpsk_symbols, norm="ortho")

# + циклический префикс
cp_ifft_symbols = np.hstack((ifft_symbols[...,-cp_len:], ifft_symbols))

TDLa_signal, impulse_res, delays_sampl = TDLa(cp_ifft_symbols.ravel(), delay_spread=delay_spread, N_fft=N, delta_f=delta_f)

noise_power = 10 ** (-signal_to_noise_ratio / 10)

# + шум
awgn_symbols = AWGN(TDLa_signal, noise_power = noise_power)

# - циклический префикс
no_cp_symbols = awgn_symbols.reshape(symbols_count, -1)[..., cp_len:]

received = np.fft.fft(no_cp_symbols, norm="ortho")

R_matr = covariance_matrix(N_fft = N, T_cp = cp_len, T_OFDM = N)
R = R_matr @ np.linalg.inv(R_matr + noise_power * np.eye(N))

h_ls = (received[0] * np.conj(pilot_symbol))
h_lmmse = R @ received[0]
W_mmse = np.conj(h_ls) / (h_ls * np.conj(h_ls) + (noise_power))
W_lmmse = np.conj(h_lmmse) / (h_lmmse * np.conj(h_lmmse) + (noise_power))

equaliz_received = received * W_mmse
equaliz_received2 = received * W_lmmse * np.exp(1j * np.pi/4)
plt.figure(1, figsize=(6,6))
plt.title(f"Delay spread = {delay_spread}, noise power = {noise_power}")
plt.scatter(np.real(received), np.imag(received), label = "received")
plt.scatter(np.real(equaliz_received), np.imag(equaliz_received), label = "equalized")
plt.scatter(np.real(equaliz_received2), np.imag(equaliz_received2), label = "equalized_lmmse")
plt.legend()
plt.grid()
plt.show()
