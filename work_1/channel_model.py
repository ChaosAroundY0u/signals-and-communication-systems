import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

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

def spatial_filter(theta_deg, phi_deg, BW_deg = 30):
    return 1.0 if ( (abs(theta_deg - BW_deg) <= BW_deg / 2) and (abs(phi_deg) <= BW_deg / 2) ) else 0.1

def transform_to_gcs(angles_deg, pointing_angles_deg):
    return angles_deg - pointing_angles_deg
    

def tworay_channel(signal):
    #cdl-a parameters
    delays_ns = np.array([0, 38.19, 40.25, 58.68])#, 46.10, 53.75, 67.08]) #nominal delay spread
    powers_db = np.array([-13.4, 0, -2.2, -4])#, -6, -8.2, -9.9])
    zod = np.array([50.2, 93.2, 93.2, 93.2])#, 94, 94, 94])
    aod = np.array([-178.1, -4.2, -4.2, -4.2])#, 90.2, 90.2, 90.2])
    #from db to linear
    powers_linear = 10 ** (powers_db / 10)
    #tdl powers (pointing direction - dominant)
    dominant_power_i = np.argmax(powers_linear)
    theta_p = zod[dominant_power_i]
    phi_p = zod[dominant_power_i]
    #applying filter
    Arxtx = np.zeros_like(zod)
    for i in range(len(zod)):
        Arxtx[i] = spatial_filter(zod[i], aod[i])
    #finally power (fr this time)
    tdl_powers = powers_linear * Arxtx
    tdl_powers /= np.max(tdl_powers)
    #impulse response
    sampling_rate = 100e6
    delays_samples = np.array([np.round(delays_ns[i] * sampling_rate / 1e9).astype(int) for i in range(len(delays_ns))])
    
    impulse_response = np.zeros(np.max(delays_samples) + 1)
    for delay, power in zip(delays_samples, powers_linear):
        if delay < len(impulse_response):
            impulse_response[delay] +=  np.sqrt(power)
            
    #for i in range(len(impulse_response)):
        
    #signal = signal / np.sqrt(1 + alpha**2)
    #delayed_signal = np.concatenate((np.zeros(delay_samples), signal[:-delay_samples])) * alpha / np.sqrt(1 + alpha**2)
    #summ = signal + delayed_signal
    signal_ans = fftconvolve(signal, impulse_response)[:5280]
    return signal_ans


# Параметры
N = 512  # количество поднесущих
symbols_count = 10
cp_len = 16
bits_count = N * symbols_count * 2  # QPSK => *2
signal_to_noise_ratio = 20

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
tworay_signal = tworay_channel(cp_ifft_symbols.ravel())


noise_power = 10 ** (-signal_to_noise_ratio / 10)

# + шум
awgn_symbols = AWGN(tworay_signal, signal_to_noise_ratio)

# - циклический префикс
no_cp_symbols = awgn_symbols.reshape(symbols_count, -1)[..., cp_len:]


received = np.fft.fft(no_cp_symbols, norm="ortho")


h_array = (received[0] / pilot_symbol)
W_mmse = np.conj(h_array) / (h_array * np.conj(h_array) + (noise_power))
#print(W_mmse)
#h_array = np.dot(np.dot(np.linalg.inv(np.dot(pilot_symbol.T, pilot_symbol)), pilot_symbol.T), received[0])
equaliz_received = received * W_mmse
#equaliz2 = received / h_array

plt.figure(2)
plt.scatter(np.real(received), np.imag(received))

plt.scatter(np.real(equaliz_received), np.imag(equaliz_received))
#plt.scatter(np.real(equaliz2), np.imag(equaliz2))
plt.grid()
plt.show()
