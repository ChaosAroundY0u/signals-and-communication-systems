clear;
close all;

%% -----------------------Variables----------------------- %%
Nfft =                                     1024;
N_with_GI =                                792; % число поднесущих с учетом Guard Interval
numSubcarrier =                            116;
symbols_count =                            5;
cp_len =                                   72;
delay_spread =                             500;
subcarrier_spacing =                       120e3;
snr =                                      70;
velocity =                                 0 * 1000/3600;
carrier_freq =                             3.5e9;
Fs =                                       Nfft * subcarrier_spacing;
%% -----------------------Main part----------------------- %%

pilot_symbol =                             bits_to_complex(1, Nfft, numSubcarrier, "QPSK");
qam_symbols =                              bits_to_complex(symbols_count-1, Nfft, numSubcarrier, "16QAM");
qam_symbols =                              [pilot_symbol, qam_symbols];

%IFFT
ifft_symbols =                             ifft(ifftshift(reshape(qam_symbols', [], symbols_count)'), [], 2);

%CP Insert
cp_ifft_symbols =                          [ifft_symbols(:, end-71:end), ifft_symbols]; % возможно стоит переделать потомушто долго прога работает

%through channel
cp_ifft_symbols =                          reshape(cp_ifft_symbols', [], 1)'; %ravel
[TDLa_signal, impulse_res, delays_sampl] = TDLa(cp_ifft_symbols, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq, symbols_count);

% + noise
awgn_symbols =                             add_noise(TDLa_signal, snr);

%CP Discard
%unravel
awgn_symbols =                             reshape(awgn_symbols', [], symbols_count)';
no_cp_symbols =                            awgn_symbols(:, 73:end);

received =                                 fftshift(fft(no_cp_symbols, [], 2));

%equalization
noise_power =                              10.^(-snr./10);
h_ls =                                     (received(1,:) .* conj(pilot_symbol));

W_mmse =                                   conj(h_ls) ./ (h_ls .* conj(h_ls) + noise_power);

equalized_received =                       received(:,:) .* W_mmse;%.* exp(1i .* pi ./ 4);

R_matr =                                   covariance_matrix(Nfft, cp_len, Nfft);
R =                                        R_matr / (R_matr + noise_power * eye(Nfft));

h_lmmse =                                  received(1, :) * R;

W_lmmse =                                  conj(h_lmmse) ./ (h_lmmse .* conj(h_lmmse) + noise_power);

equalized_received2 =                      received(:, :) .* W_lmmse .* exp(1i .* pi ./ 4);

%Plots
figure;
%scatter(real(received), imag(received), 'red');
scatter(real(equalized_received(:, 117:792+116)), imag(equalized_received(:, 117:792+116)), 'blue');
hold on;
scatter(real(qam_symbols(116+1:792+116)), imag(qam_symbols(116+1:792+116)),'red');

scatter(real(equalized_received2(116+1:792+116)), imag(equalized_received2(116+1:792+116)), 'green');

grid on;
%legend('received', 'equlized-ls', 'equalized-lmmse');
%figure; plot_spectrum2(cp_ifft_symbols(1, :), Fs, "red", length(cp_ifft_symbols(1, :)));

%% -----------------------Functions----------------------- %%
function plot_spectrum2(signal, fs, colour, N)
    grid on;
    f = (-N/2:N/2-1) * (fs / N);
    spectrum = fftshift(fft(signal, [], 2));
    plot(f, 20*log10(spectrum), color = colour);
    xlabel('Frequency, Hz');
    ylabel('Magnitude, dB');
    xlim([-fs/2, fs/2]);
    hold on;
end

function bpsk_out = bpsk_modulate(x)
    bpsk_out = 1./sqrt(2) .* ((1-2*x(1:end)) + 1j .* (1-2*x(1:end)));
end

function qpsk_out = qpsk_modulate(x)
    qpsk_out =  1./sqrt(2) .* (1-2*x(1:2:end) + 1j * (1-2*x(2:2:end)));
end

function qam16_out = qam16_modulate(x)
    qam16_out = 1./sqrt(10) .* ( (1-2.*x(1:4:end)) .* (2-(1-2.*x(3:4:end))) + 1j .* (1-2.*x(2:4:end)) .* (2-(1-2*x(4:4:end))) );
end

function qam64_out = qam64_modulate(x)
    qam64_out = 1/sqrt(42) .* ( (1-2*x(1:6:end)) .* ( 4 - (1 - 2*x(3:6:end)).*(2 - (1 - 2*x(5:6:end))) ) + 1j * (1 - 2*x(2:6:end)) .* (4 - (1 - 2*x(4:6:end)) .* (2 - (1 - 2*x(6:6:end)))) );
end

function out = bits_to_complex(symbs_count, nfft, num_subcarrier, modulation_technique)
    modulations =                                   ["BPSK", "QPSK", "", "16QAM","", "64QAM"];
    scale_var =                                     find(modulations == modulation_technique);
    x_count =                                       nfft * symbs_count * scale_var; % число бит
    x =                                             randi([0, 1], 1, x_count);
    out =                                           zeros(1, symbs_count .* nfft);
    switch scale_var
        case 1
            modulate =                              @bpsk_modulate;
        case 2
            modulate =                              @qpsk_modulate;
        case 4
            modulate =                              @qam16_modulate;
        case 6
            modulate =                              @qam64_modulate;
    end
    for iterable = 1:symbs_count
        x_use =                                     x(iterable*scale_var*nfft - scale_var*nfft + scale_var*num_subcarrier + 1:iterable*scale_var*nfft - scale_var * num_subcarrier);
        %symbol =                                    x(iterable*scale_var*nfft-scale_var*nfft+1:iterable*scale_var*nfft);
        modulated =                                 modulate(x_use);
        out(iterable*nfft-nfft+1:iterable*nfft) =   [zeros(1, num_subcarrier), modulated, zeros(1, num_subcarrier)];
    end
end

function noise_added = add_noise(signal, SNR)
    noise_power = 10 .^ (-SNR./10);
    noise = sqrt(noise_power ./ 2) .* (randn(size(signal)) + 1j .* randn(size(signal)));
    noise_added = signal + noise;
end

function [output_signal, powers_linear, delays_samples] = TDLa(signal, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq, num_symbs)
    delays_ns =             [0 .3819 .4025 .5868 .4610 .5375 .6708] .* delay_spread;
    powers_db =             [-13 0 -2.2 -4 -6 -8.2 -9.9];
    powers_linear =         10 .^ (powers_db ./ 10);
    %powers_linear =         powers_linear ./ sum(powers_linear);
    sample_rate =           Nfft .* subcarrier_spacing;
    Ts =                    1 / sample_rate; % symbol time
    delays_samples =        cast(round(delays_ns .* sample_rate / 1e9), 'int32');
    
    f_d =                   (velocity .* carrier_freq) ./ (3*1e8);
    fs_channel =            max(10 * f_d, 1000);

    num_samples =           length(signal);
    num_taps =              length(delays_samples);
    time_process =          num_samples / sample_rate; % время прохождения сигнала через канал
    num_channel_samples =   ceil(time_process * fs_channel);
    t_channel =             (0:num_channel_samples-1) / fs_channel;
    
    N0 =                    5; % число синусоид
    beta =                  pi * (1:N0) / N0;
    alpha =                 (1:N0) * (2*pi/N0);
    
    channel_taps_history =  zeros(num_taps, num_channel_samples);
    
    for x = 1:num_taps
        in_phase =  zeros(1, num_channel_samples);
        quadrture = zeros(1, num_channel_samples);
        
        for n = 1:N0
            omega_n =     2*pi*f_d * cos(alpha(n));
            phase_shift = 2*pi*rand();  % случайная фаза для каждого тапа
            in_phase =    in_phase + cos(beta(n) + phase_shift) .* cos(omega_n * t_channel);
            quadrture =   quadrture + sin(beta(n) + phase_shift) .* cos(omega_n * t_channel);
        end
        
        fading =    (in_phase + 1j * quadrture) / sqrt(N0);
        
        % Базовый коэффициент тапа
        base_tap =                   sqrt(powers_linear(x)/2) * (randn(1,1) + 1j * randn(1,1));
        channel_taps_history(x, :) = base_tap .* fading;
    end
    
    max_lag = 70; % samples
    autocorr(max_lag, channel_taps_history(1, :), time_process, num_symbs, sample_rate, fs_channel, f_d, velocity, N0);
    
    output_signal = zeros(1, num_samples);
    max_delay = max(delays_samples);
    
    samples_per_channel_sample = round(sample_rate / fs_channel);
    
    fprintf('Информация о дискретизации:\n');
    fprintf('  Частота дискретизации сигнала: %.2f Гц\n', sample_rate);
    fprintf('  Частота обновления канала: %.2f Гц\n', fs_channel);
    fprintf('  1 отсчет канала используется для %.0f отсчетов сигнала\n', samples_per_channel_sample);
    
    for n_signal = 1:num_samples
        % определяем, какой отсчет канала использовать
        % время текущего отсчета сигнала:
        t_signal = (n_signal - 1) / sample_rate;
        
        % индекс отсчета канала для этого времени
        n_channel = floor(t_signal * fs_channel) + 1;
        n_channel = min(n_channel, num_channel_samples); % защита от выхода за границы
        
        current_output = 0;
        
        % 2. Свертка с текущим каналом
        for tap_idx = 1:num_taps
            delay = delays_samples(tap_idx);
            
            if n_signal > delay  % проверка, что не выходим за границы
                % берём значение тапа для текущего момента времени
                tap_value = channel_taps_history(tap_idx, n_channel);
                
                % берём отсчет сигнала с соответствующей задержкой
                signal_value = signal(n_signal - delay);
                
                current_output = current_output + tap_value * signal_value;
            end
        end
        
        output_signal(n_signal) = current_output;
    end
end

function plot_autocorr = autocorr(max_lag, process, time_process, num_symbs, sample_rate, fs_channel, f_d, velocity, N0)
    lags = 0:max_lag-1;
    corr = zeros(1, max_lag);
    N = length(process);
    for k = 1:max_lag
        idx = k - 1;
        corr(k) = mean(process(1:N-idx) .* conj(process(1 + idx:N)));
    end
    time_lags = lags ./ sample_rate;
    figure;
    corr = corr ./ max(corr);
    tau = (0:max_lag-1) ./ fs_channel;
    J = besselj(0, 2*pi*f_d*tau);
    %J = J ./ max(J);
    plot(time_lags .*1e5, corr);
    hold on;
    plot(time_lags .*1e5, J);
    xlabel('Сдвиг, мкс');
    ylabel('корр(сдвиг)');
    T = get_correlation_time_simple(corr, J, sample_rate) *1e5;
    title(sprintf('время процесса %.6f c (%.i символов)', time_process, num_symbs), sprintf('скорость приемника = %.2f км/ч, число синусоид = %d, время корреляции = %.2f мкс', velocity .*3600./1000, N0, T) );
    legend('вычисляемая акф', 'ф-я Бесселя 0 порядка')
    grid on;
end

function T_corr = get_correlation_time_simple(signal1, signal2, fs)
    
    x1 = signal1 - mean(signal1);
    x2 = signal2 - mean(signal2);
    [corr_vals, ~] = xcorr(x1, x2, 'normalized');
    
    center_idx = ceil(length(corr_vals)/2);
    corr_peak = corr_vals(center_idx);
    
    threshold = abs(corr_peak) * 0.2;
    right_side = abs(corr_vals(center_idx:end));
    
    % находим первый индекс ниже порога
    idx = find(right_side < threshold, 1);
    
    if isempty(idx)
        T_corr = length(signal1)/fs;  % если не упало ниже порога
    else
        T_corr = (idx - 1) / fs;
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

