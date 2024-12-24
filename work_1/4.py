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
    target_noise_watts = 10 ** (signal_to_noise_ratio / 10)
    signal_power = np.mean(np.abs(signal) ** 2)
    noise_power = signal_power / target_noise_watts
    noise = np.sqrt(noise_power) * (np.random.normal(size=signal.shape) + 1j * np.random.normal(size=signal.shape))
    return signal + noise

def tworay_channel(signal, delay_samples, alpha):
    signal = signal / np.sqrt(1 + alpha**2)
    delayed_signal = np.concatenate((np.zeros(delay_samples), signal[:-delay_samples])) * alpha / np.sqrt(1 + alpha**2)
    summ = signal + delayed_signal
    return summ

def split_array(arr, n):
    # Вычисляем размер каждой части
    k, m = divmod(len(arr), n)
    return [arr[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


# Параметры
N = 512  # количество поднесущих
symbols_count = 10
cp_len = 16
bits_count = N * symbols_count * 2  # QPSK => *2
signal_to_noise_ratio = 50
data_bits = np.random.randint(0, 2, bits_count)

# QPSK
qpsk_symbols = qpsk_modulate(data_bits)

pilot_symbol = qpsk_symbols[:N] # первый символ из 10 - пилотный
qpsk_symbols = qpsk_symbols[N:] # 9 символов после первого

# IFFT
ifft_symbols = np.fft.ifft(qpsk_symbols)

# + циклический префикс
cp_ifft_symbols = np.hstack((ifft_symbols[-cp_len:], ifft_symbols))

# + шум
awgn_symbols = AWGN(cp_ifft_symbols, signal_to_noise_ratio)

# - циклический префикс
no_cp_symbols = awgn_symbols[cp_len:]

delay_sampl = 10
tworay_signal = tworay_channel(no_cp_symbols, delay_sampl, 0.5)

received = np.fft.fft(tworay_signal)

#тут разбиваем ресивд шоб оценку посчитать (делим на пилота)
received_chunks = split_array(received, symbols_count - 1)
#теперь делим 9 полученых чанков на пилот
h_array = []
for i in range(len(received_chunks)):
    h = np.array(received_chunks[i]) / np.array(pilot_symbol)
    h_array.append(h)
    
#делим ресивд чанки на оценку (эквализация)
equalization_arr = []
for i in range(len(received_chunks)):
    equalization = np.array(received_chunks[i]) / h_array[i]
    equalization_arr.append(equalization)

#тут я двойной массив в один превращаю
equaliz_received = [item for sublist in equalization_arr for item in sublist]
print(equaliz_received)

plt.figure(2)
plt.scatter(np.real(received), np.imag(received))
plt.grid()
plt.show()

plt.figure(3)
plt.scatter(np.real(equaliz_received), np.imag(equaliz_received))
plt.grid()
plt.show()
