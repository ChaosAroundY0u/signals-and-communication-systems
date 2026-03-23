function noise_added = add_noise(signal, SNR)
    noise_power = 10 .^ (-SNR./10);
    noise = sqrt(noise_power ./ 2) .* (randn(size(signal)) + 1j .* randn(size(signal)));
    noise_added = signal + noise;
end
