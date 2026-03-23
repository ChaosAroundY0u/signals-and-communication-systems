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
