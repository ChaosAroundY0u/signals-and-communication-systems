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
    T = correlation_time(corr, J, sample_rate) *1e5;
    %title(sprintf('время процесса %.6f c (%.i символов)', time_process, num_symbs), sprintf('скорость приемника = %.2f км/ч, число синусоид = %d, время корреляции = %.2f мкс', velocity .*3600./1000, N0, T) );
    legend('вычисляемая акф', 'ф-я Бесселя 0 порядка')
    grid on;
end
