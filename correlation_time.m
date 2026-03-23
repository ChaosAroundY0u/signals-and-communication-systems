function T_corr = correlation_time(signal1, signal2, fs)
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
