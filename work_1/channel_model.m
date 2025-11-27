clear;
close all;

%% -----------------------Variables----------------------- %%
Nfft = 1024;
symbols_count = 2;
cp_len = 72;
delay_spread = 900;
subcarrier_spacing = 15e3;
bits_count = Nfft * symbols_count * 2;
snr = 30;
bits = randi([0, 1], 1, bits_count);
velocity = 0.1 .* 3*10e8 * 1000/3600;
carrier_freq = 3.5e9;
%% -----------------------Main part----------------------- %%
%qpsk modulate
qpsk_symbols = zeros(1, symbols_count .* Nfft);
for i = 1:symbols_count
    symbol = bits(i*2*Nfft-2*Nfft+1:i*2*Nfft);
    modulated = qpsk_modulate(symbol);
    qpsk_symbols([i.*Nfft-Nfft+1:i.*Nfft]) = modulated;
end
%pilot symbol
pilot_symbol = qpsk_symbols(1:Nfft);

%IFFT
ifft_symbols = ifft(ifftshift(reshape(qpsk_symbols', [], symbols_count)'), [], 2);

%CP Insert
cp_ifft_symbols = [ifft_symbols(:, end-71:end), ifft_symbols]; % возможно стоит переделать потомушто долго прога работает

%through channel
cp_ifft_symbols = reshape(cp_ifft_symbols', [], 1)'; %ravel
[TDLa_signal, impulse_res, delays_sampl] = TDLa(cp_ifft_symbols, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq);
% 
% % + noise
% awgn_symbols = add_noise(TDLa_signal, snr);
% 
% %CP Discard
% %unravel
% awgn_symbols = reshape(awgn_symbols', [], symbols_count)';
% no_cp_symbols = awgn_symbols(:, 73:end);
% 
% received = fftshift(fft(no_cp_symbols, [], 2));
% 
% %equalization
% noise_power = 10.^(-snr./10);
% h_ls = (received(1, :) .* conj(pilot_symbol));
% W_mmse = conj(h_ls) ./ (h_ls .* conj(h_ls) + noise_power);
% equalized_received = received .* W_mmse;
% 
% R_matr = covariance_matrix(Nfft, cp_len, Nfft);
% R = R_matr / (R_matr + noise_power * eye(Nfft));
% 
% h_lmmse = received(1, :) * R;
% W_lmmse = conj(h_lmmse) ./ (h_lmmse .* conj(h_lmmse) + noise_power);
% equalized_received2 = received .* W_lmmse .* exp(1i .* pi ./ 4);
% 
% 
% %Plots
% figure;
% scatter(real(received), imag(received), 'red');
% hold on;
% scatter(real(equalized_received), imag(equalized_received), 'blue');
% scatter(real(equalized_received2), imag(equalized_received2), 'green');
% grid on;
% legend('received', 'equlized', 'equalized-lmmse');


%% -----------------------Functions----------------------- %%
function qpsk_out = qpsk_modulate(x)
    %return 1/np.sqrt(2) * ((1-2*data_bits[0::2]) + 1j*(1-2*data_bits[1::2])) #TS 38.211 - 5.1.3
    qpsk_out =  1./sqrt(2) .* (1-2*x(1:2:end) + 1j * (1-2*x(2:2:end)));
end

function noise_added = add_noise(signal, SNR)
    noise_added = awgn(signal, SNR);
end

function [output_signal, powers_linear, delays_samples] = TDLa(signal, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq)
    delays_ns = [0 .3819 .4025 .5868 .4610 .5375 .6708] .* delay_spread;
    powers_db = [-13 0 -2.2 -4 -6 -8.2 -9.9];
  
    powers_linear = 10 .^ (powers_db ./ 10);
    powers_linear = powers_linear ./ sum(powers_linear);
    sample_rate = Nfft .* subcarrier_spacing;
    Ts = 1 / sample_rate; % период дискретизации
    delays_samples = cast(round(delays_ns .* sample_rate / 1e9), 'int32');
        
    num_samples = length(signal);
    num_taps = length(delays_samples);
    
    % Доплеровская частота
    c = 3e8; % скорость света
    % New model (N = 12)
    N = [1 2 3 4 5 6 7 8 9 10 11 12];
    N0 = length(N) / 4;
    alpha_n = N.*(2.*pi./length(N));
    beta_n = N(1:N0).*pi./N0;
    omega_m = 2.*pi.*carrier_freq.*velocity./c; % maximum Doppler shift;
    omega_n = omega_m .* cos(alpha_n);
    
    channel_taps_history = zeros(num_taps, num_samples);
    
    t = (0:num_samples-1) * Ts;
    
    for x = 1:num_taps
        waveform = 0;
        T_t = zeros(1, length(t));
        for g = 1:length(t)
            waveform = 0;
            theta_n = 2.*pi.*rand(1, 3);
            for i = 1:N0
                waveform = waveform + sqrt(2./N0) .* (cos(beta_n(i)) + 1j.*sin(beta_n(i)));%.*cos(omega_n(i) + theta_n(i));
            end
            T_t(g) = waveform .* cos(sum(omega_n).*t(g) + sum(theta_n));       
        end
        
        % Базовый коэффициент тапа + доплеровская модуляция
        base_tap = sqrt(powers_linear(i)) * (randn(1, 1) + 1j * randn(1, 1)) / sqrt(2);
        channel_taps_history(x, :) = base_tap .* T_t;
    end
    
    % Свёртка с изменяющимся каналом
    output_signal = zeros(1, num_samples);
    max_delay = max(delays_samples);
    
    for n = 1:num_samples
        % Импульсная характеристика для текущего момента времени
        impulse_response = zeros(1, max_delay + 1);
        
        for i = 1:num_taps
            delay = delays_samples(i);
            if delay + 1 <= length(impulse_response)
                impulse_response(delay + 1) = impulse_response(delay + 1) + channel_taps_history(i, n);
            end
        end
        
        % Свёртка для текущего отсчёта
        if n <= length(impulse_response)
            segment = signal(1:n);
            ir_segment = impulse_response(1:n);
            output_signal(n) = sum(segment .* fliplr(ir_segment));
        else
            segment_start = n - length(impulse_response) + 1;
            segment_end = n;
            segment = signal(segment_start:segment_end);
            output_signal(n) = sum(segment .* fliplr(impulse_response));
        end
    end
end



function R_hh = covariance_matrix(Nfft, T_cp, T_OFDM)
    alpha = T_cp ./ T_OFDM;
    R_hh = zeros(Nfft);
    for k = 1:Nfft
        for m = 1:Nfft
            if k == m
                R_hh(k, m) = 1;
            else
                delta = k - m;
                phase = -pi .* delta .* alpha;
                R_hh(k, m) = (1 / alpha) .* sinc(delta .* alpha) .* exp(1i .* phase);

            end
        end
    end
end
