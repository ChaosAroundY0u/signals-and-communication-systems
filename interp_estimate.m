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
