function [output_signal, powers_linear, delays_samples] = TDLa(signal, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq, num_symbs)
    delays_ns =             [0 .3819 .4025 .5868 .4610 .5375 .6708] .* delay_spread;
    powers_db =             [-13 0 -2.2 -4 -6 -8.2 -9.9];
    powers_linear =         10 .^ (powers_db ./ 10);
    %powers_linear =         powers_linear ./ sum(powers_linear);
    sample_rate =           Nfft .* subcarrier_spacing;
    delays_samples =        cast(round(delays_ns .* sample_rate / 1e9), 'int32');

    % ЧАСТОТА ДИСКРЕТИЗАЦИИ КАНАЛА должна быть >= частоты Доплера
    f_d =                   (velocity .* carrier_freq) ./ (3*1e8);
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
