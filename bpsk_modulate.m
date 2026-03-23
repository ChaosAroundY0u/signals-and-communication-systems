function bpsk_out = bpsk_modulate(x)
    bpsk_out = 1./sqrt(2) .* ((1-2*x(1:end)) + 1j .* (1-2*x(1:end)));
end
