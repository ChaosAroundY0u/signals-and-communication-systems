clear;
close all;

%% -----------------------Variables----------------------- %%
Nfft =                                     1024;
N_with_GI =                                792; % число поднесущих с учетом Guard Interval
numSubcarrier =                            116;
symbols_count =                            12;
cp_len =                                   72;
delay_spread =                             500;
subcarrier_spacing =                       120e3;
snr =                                      70;
velocity =                                 60 * 1000/3600;
carrier_freq =                             3.5e9;
Fs =                                       Nfft * subcarrier_spacing;
%% -----------------------Main part----------------------- %%

%pilot_symbol =                             bits_to_complex(1, Nfft, numSubcarrier, 2);
symb_counter =   0;
pilots_counter = 0;
symbols =       [];
pilots =        [];
while symb_counter < symbols_count
    symbol = bits_to_complex(1, Nfft, numSubcarrier, 4);
    if mod(symb_counter, 10) == 0
        symbol = bits_to_complex(1, Nfft, numSubcarrier, 2);
        pilots = [pilots, symbol];
        pilots_counter = pilots_counter + 1;
    end
    symbols = [symbols, symbol];
    symb_counter = symb_counter + 1;
end

pilots =                                   reshape(pilots', [], pilots_counter)';
qam_symbols =                              symbols;
% IFFT
ifft_symbols =                             ifft(ifftshift(reshape(qam_symbols', [], symbols_count)'), [], 2);
%CP Insert
cp_ifft_symbols =                          [ifft_symbols(:, end-71:end), ifft_symbols];
%Through channel
cp_ifft_symbols =                          reshape(cp_ifft_symbols', [], 1)'; %ravel
[TDLa_signal, impulse_res, delays_sampl] = TDLa(cp_ifft_symbols, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq, symbols_count);
%+ noise
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

W_mmse =                                   conj(h_ls_estimates) ./ (h_ls_estimates .* conj(h_ls_estimates) + noise_power);
equalized_received = [];
for i = 1:symbols_count
    estimate = W_mmse(fix((i-1) / 10) + 1, :);
    equalized = received(i, :) .* estimate;
    equalized_received = [equalized_received; equalized];
end
R_matr =                                   covariance_matrix(Nfft, cp_len, Nfft);
R =                                        R_matr / (R_matr + noise_power * eye(Nfft));

h_lmmse =                                  h_ls_estimates * R;
W_lmmse =                                  conj(h_lmmse) ./ (h_lmmse .* conj(h_lmmse) + noise_power);

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

figure;
scatter(real(received(:, 117:792+116)), imag(received(:, 117:792+116)), "r"); hold on;
scatter(real(equalized_received(:, 117:792+116)), imag(equalized_received(:, 117:792+116)), "b"); hold on;
scatter(real(equalized_received2(:, 116+1:792+116)), imag(equalized_received2(:, 116+1:792+116)), "g"); hold on;
scatter(real(equalized_received3(:, 116+1:792+116)), imag(equalized_received3(:, 116+1:792+116)), "m"); grid on;
