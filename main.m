clear;
close all;


%% -----------------------Variables----------------------- %%
Nfft =                                     64;
N_with_GI =                                64; % число поднесущих с учетом Guard Interval
numSubcarrier =                            0;  % число поднесущих из GI
symbols_count =                            100;
cp_len =                                   16;
delay_spread =                             200;
subcarrier_spacing =                       312500;
snr =                                      30;
velocity =                                 10 * 1000/3600;
carrier_freq =                             2.4e9;
Fs =                                       Nfft * subcarrier_spacing;
%% -----------------------Main part----------------------- %%

symb_counter =   0;
pilots_counter = 0;
symbols =       [];
pilots =        [];
[symbols, encoded] = bits_to_complex(symbols_count, Nfft, numSubcarrier, 2);

symb_list =                                [1:symbols_count];
idx =                                      (mod(symb_list, 10) == 1) .*symb_list;
idx =                                      idx(~idx == 0);
symbols =                                  reshape(symbols', [], symbols_count)';
pilots =                                   symbols(idx, :);
qam_symbols =                              symbols;
% IFFT
ifft_symbols =                             ifft(ifftshift(reshape(qam_symbols', [], symbols_count)'), [], 2);
%CP Insert
cp_len = 2.*cp_len;
cp_ifft_symbols =                          [ifft_symbols(:, end-(cp_len-1):end), ifft_symbols];
%Through channel
cp_ifft_symbols =                          reshape(cp_ifft_symbols', [], 1)'; %ravel
[TDLa_signal, impulse_res, delays_sampl] = TDLa(cp_ifft_symbols, delay_spread, Nfft, subcarrier_spacing, velocity, carrier_freq, symbols_count);
%+ noise
awgn_symbols =                             add_noise(TDLa_signal, snr);
%CP Discard
%unravel
awgn_symbols =                             reshape(awgn_symbols', [], symbols_count)';
no_cp_symbols =                            awgn_symbols(:, (cp_len+1):end);
received =                                 fftshift(fft(no_cp_symbols, [], 2));
%equalization
noise_power =                              10.^(-snr./10);

h_ls_estimates = [];
for i = 1:length(idx)
    h_ls = received(idx(i), :) .* conj(pilots(i, :));
    h_ls_estimates = [h_ls_estimates, h_ls];
end
   h_ls_estimates = reshape(h_ls_estimates', [], length(idx))';

W_mmse =                                   conj(h_ls_estimates) ./ (h_ls_estimates .* conj(h_ls_estimates) + noise_power);
equalized_received = [];
for i = 1:symbols_count
    estimate = W_mmse(fix((i-1) / 10) + 1, :);
    equalized = received(i, :) .* estimate;
    equalized_received = [equalized_received; equalized];
end
R_matr =                                   covariance_matrix(Nfft, cp_len, Nfft);
R =                                        R_matr / (R_matr + noise_power * eye(2*Nfft));

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
symb_list = [1:symbols_count];
idx = ~(mod(symb_list, 10) == 1) .* symb_list;
idx = idx(~idx == 0);
scatter(real(received(idx, 1:64)), imag(received(idx, 1:64)), "r"); hold on;
scatter(real(equalized_received(idx, 1:64)), imag(equalized_received(idx, 1:64)), "b"); hold on;
scatter(real(equalized_received2(idx, 1:64)), imag(equalized_received2(idx, 1:64)), "g"); hold on;
scatter(real(equalized_received3(idx, 1:64)), imag(equalized_received3(idx, 1:64)), "m"); grid on;

demapping = hard_demapping(equalized_received3);
