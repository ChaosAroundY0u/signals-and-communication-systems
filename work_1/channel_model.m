
%% -----------------------Variables----------------------- %%
Nfft =                                     1024;
N_with_GI =                                1024; % число поднесущих с учетом Guard Interval
numSubcarrier =                            0;
symbols_count =                            50;
cp_len =                                   72;
delay_spread =                             500;
subcarrier_spacing =                       120e3;
snr =                                      60;
velocity =                                 60 * 1000/3600;
carrier_freq =                             3.5e9;
Fs =                                       Nfft * subcarrier_spacing;
%коэффициенты полиномов идут слева направо, т.е.:
g1 =                                       [1 0]; % 0*x + 1 или  вернее 1 + 0*x
g2 =                                       [1 1]; % 1*x + 1 -/-/-/-/-/- 1 + 1*x
%под знаком сложения подразумевается операция xor (сложение по модулю 2)
%% -----------------------Main part----------------------- %%

symb_counter =   0;
pilots_counter = 0;
symbols =       [];
pilots =        [];
while symb_counter < symbols_count
    symbol = bits_to_complex(1, Nfft, numSubcarrier, "16QAM");
    if mod(symb_counter, 10) == 0
        symbol = bits_to_complex(1, Nfft, numSubcarrier, "QPSK");
        pilots = [pilots, symbol];
        pilots_counter = pilots_counter + 1;
    end
    symbols = [symbols, symbol];
    symb_counter = symb_counter + 1;
end

pilots = reshape(pilots', [], pilots_counter)';

%qpsk_symbols =                             bits_to_complex(symbols_count, Nfft, numSubcarrier, "QPSK");
%pilot_symbol =                             bits_to_complex(1, Nfft, numSubcarrier, "QPSK", g1, g2);
%pilot_symbol =                             bits_to_complex(1, Nfft, numSubcarrier, "QPSK");
%qam_symbols =                              bits_to_complex(symbols_count-1, Nfft, numSubcarrier, "QPSK", g1, g2);
%qam_symbols =                              bits_to_complex(symbols_count-1, Nfft, numSubcarrier, "16QAM");
%qam_symbols =                              [pilot_symbol, qam_symbols];
qam_symbols =                               symbols;

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

h_ls_estimates = [];
for i = 1:pilots_counter
    h_ls = received((i-1).*10+1, :) .* conj(pilots(i, :));
    h_ls_estimates = [h_ls_estimates, h_ls];
end
h_ls_estimates = reshape(h_ls_estimates', [], pilots_counter)';
h_ls =                                     (received(1,:) .* conj(pilots(1,:))); %conj(pilot_symbol));

W_mmse =                                   conj(h_ls_estimates) ./ (h_ls_estimates .* conj(h_ls_estimates) + noise_power);
equalized_received = [];
for i = 1:symbols_count
    estimate = W_mmse(fix((i-1) / 10) + 1, :);
    equalized = received(i, :) .* estimate;
    equalized_received = [equalized_received; equalized];
end

%equalized_received =                       received(:,:) .* W_mmse;%.* exp(1i .* pi ./ 4);

R_matr =                                   covariance_matrix(Nfft, cp_len, Nfft);
R =                                        R_matr / (R_matr + noise_power * eye(Nfft));

h_lmmse =                                  h_ls_estimates * R;
W_lmmse =                                  conj(h_lmmse) ./ (h_lmmse .* conj(h_lmmse) + noise_power);

%equalized_received2 =                      received(:, :) .* W_lmmse;
equalized_received2 = [];
for i = 1:symbols_count
    estimate = W_lmmse(fix((i-1) / 10) + 1, :);
    equalized = received(i, :) .* estimate;
    equalized_received2 = [equalized_received2; equalized];
end

H_interp = interp_estimate(h_ls_estimates, symbols_count, snr, Nfft, subcarrier_spacing, velocity, carrier_freq);

equalized_received3 = [];
for i = 1:symbols_count
    equalized3 = received(i, :) ./ H_interp(i, :);
    equalized_received3 = [equalized_received3; equalized3];
end
%Plots
figure;
%scatter(real(received(:, 117:792+116)), imag(received(:, 117:792+116)), 'red');
hold on;
scatter(real(equalized_received(:, 117:792+116)), imag(equalized_received(:, 117:792+116)), 'blue');
hold on;
%scatter(real(qam_symbols(116+1:792+116)), imag(qam_symbols(116+1:792+116)),'red');

scatter(real(equalized_received2(:, 116+1:792+116)), imag(equalized_received2(:, 116+1:792+116)), 'green');
hold on;
scatter(real(equalized_received3(:, 116+1:792+116)), imag(equalized_received3(:, 116+1:792+116)), 'magenta');

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
    %coded_bits_per_symbol =                         nfft*scale_var*2+2;
    %[x, coded_bits] =                               conv_encoder(x, g1, g2);
    %x =                                             [x, zeros(1, symbs_count * coded_bits_per_symbol - coded_bits)];
    %out =                                           [zeros(1, (symbs_count .* nfft) * 2 + 1*(symbs_count-1) + 1), zeros(1, symbs_count * coded_bits_per_symbol - coded_bits)];
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
        %x_use =                                     x(((iterable*scale_var*nfft) - (scale_var*nfft + scale_var*num_subcarrier))*2+ 2*(iterable-1) + 1:(iterable*scale_var*nfft - scale_var * num_subcarrier)*2+2*iterable);
        x_use =                                     x(iterable*scale_var*nfft - scale_var*nfft + scale_var*num_subcarrier + 1:iterable*scale_var*nfft - scale_var * num_subcarrier);
        %symbol =                                    x(iterable*scale_var*nfft-scale_var*nfft+1:iterable*scale_var*nfft);
        modulated =                                 modulate(x_use);
        %out((iterable*nfft-nfft)*2 + 1*(iterable-1)+1:iterable*nfft*2 + 1*(iterable-1) + 1) =   [zeros(1, num_subcarrier), modulated, zeros(1, num_subcarrier)];
        out(iterable*nfft-nfft+1:iterable*nfft) =   [zeros(1, num_subcarrier), modulated, zeros(1, num_subcarrier)];
    end
end

function noise_added = add_noise(signal, SNR)
    noise_power = 10 .^ (-SNR./10);
    %noise = np.sqrt(noise_power / 2) * (np.random.normal(size=signal.shape) + 1j * np.random.normal(size=signal.shape))
    %return signal + noise
    noise = sqrt(noise_power ./ 2) .* (randn(size(signal)) + 1j .* randn(size(signal)));
    %noise_added = awgn(signal, SNR);
    noise_added = signal + noise;
end

function [output_signal, powers_linear, delays_samples] = TDLa(signal, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq, num_symbs)
    delays_ns =             [0 .3819 .4025 .5868 .4610 .5375 .6708] .* delay_spread;
    powers_db =             [-13 0 -2.2 -4 -6 -8.2 -9.9];
    powers_linear =         10 .^ (powers_db ./ 10);
    %powers_linear =         powers_linear ./ sum(powers_linear);
    sample_rate =           Nfft .* subcarrier_spacing;
    delays_samples =        cast(round(delays_ns .* sample_rate / 1e9), 'int32');
    
    % ЧАСТОТА ДИСКРЕТИЗАЦИИ КАНАЛА должна быть >= частоты Доплера
    f_d =                   (velocity .* carrier_freq) ./ (3*1e8);  % исправлено: 3e8, а не 3*10e8
    fs_channel =            max(10 * f_d, 1); 

    num_samples =           length(signal);
    num_taps =              length(delays_samples);
    time_process =          num_samples / sample_rate; % время прохождения сигнала через канал
    num_channel_samples =   ceil(time_process * fs_channel);
    t_channel =             (0:num_channel_samples-1) / fs_channel;
    
    % Генерация доплеровских замираний (модель Джейкса)
    N0 =                    5; % число синусоид
    beta =                  pi * (1:N0) / N0;
    alpha =                 (1:N0) * (2*pi/N0);
    
    channel_taps_history =  zeros(num_taps, num_channel_samples);
    
    % Генерация ОДНОГО процесса замираний для всех тапов (но с разными фазами)
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
    
    % При v=0 убедимся, что канал постоянен
    % if velocity == 0
    %     for x = 1:num_taps
    %         channel_taps_history(x, :) = channel_taps_history(x, 1);  % все отсчеты одинаковые
    %     end
    % end
    max_lag = 70; % samples
    autocorr(max_lag, channel_taps_history(1, :), time_process, num_symbs, sample_rate, fs_channel, f_d, velocity, N0);
    
    % ====== ИСПРАВЛЕННАЯ СВЕРТКА ======
    output_signal = zeros(1, num_samples);
    max_delay = max(delays_samples);
    
    samples_per_channel_sample = round(sample_rate / fs_channel);
    
    fprintf('Информация о дискретизации:\n');
    fprintf('  Частота дискретизации сигнала: %.2f Гц\n', sample_rate);
    fprintf('  Частота обновления канала: %.2f Гц\n', fs_channel);
    fprintf('  1 отсчет канала используется для %.0f отсчетов сигнала\n', samples_per_channel_sample);
    
    for n_signal = 1:num_samples
        % 1. Определяем, какой отсчет канала использовать
        % Время текущего отсчета сигнала:
        t_signal = (n_signal - 1) / sample_rate;
        
        % Индекс отсчета канала для этого времени:
        n_channel = floor(t_signal * fs_channel) + 1;
        n_channel = min(n_channel, num_channel_samples); % защита от выхода за границы
        
        current_output = 0;
        
        % 2. Свертка с текущим каналом
        for tap_idx = 1:num_taps
            delay = delays_samples(tap_idx);
            
            if n_signal > delay  % проверка, что не выходим за границы
                % Берём значение тапа для текущего момента времени
                tap_value = channel_taps_history(tap_idx, n_channel);
                
                % Берём отсчет сигнала с соответствующей задержкой
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
    % Нормализованная кросс-корреляция
    x1 = signal1 - mean(signal1);
    x2 = signal2 - mean(signal2);
    [corr_vals, ~] = xcorr(x1, x2, 'normalized');
    
    center_idx = ceil(length(corr_vals)/2);
    corr_peak = corr_vals(center_idx);
    
    threshold = abs(corr_peak) * 0.2;
    right_side = abs(corr_vals(center_idx:end));
    
    % Находим первый индекс ниже порога
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



function H_interp = interp_estimate(estimates, num_symbs, snr, nfft, scs, velocity, carrier_freq)

    pilot_indices = [1:10:num_symbs];
    target_indices = [1:num_symbs];
    snr_linear = 10^(snr/10);
    sigma_square = 1 / snr_linear;
    fd = (velocity .* carrier_freq) ./ (3*1e8);
    sample_rate = nfft .* scs;
    symb_time = nfft ./ sample_rate;

    time_all = (0:num_symbs-1) * symb_time;
    time_pilots = time_all(pilot_indices);
    time_targets = time_all(target_indices);
    
    tau_pp = abs(time_pilots' - time_pilots);
    R_pp = besselj(0, 2 * pi * fd * tau_pp);
    R_pp = R_pp + sigma_square * eye(length(pilot_indices));

    tau_tp = abs(time_targets' - time_pilots);
    R_tp = besselj(0, 2 * pi * fd * tau_tp);
    weights = R_tp / R_pp;

    H_interp = [];
    H_sum = 0;
    for i = 1:num_symbs
        for j = 1:size(weights, 2)
            H = weights(i, j) .* estimates(j, :);
            H_sum = H_sum + H;
        end
        H_interp = [H_interp; H_sum];
        H_sum = 0;
    end
end
