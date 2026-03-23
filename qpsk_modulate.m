function qpsk_out = qpsk_modulate(x)
    qpsk_out =  1./sqrt(2) .* (1-2*x(1:2:end) + 1j * (1-2*x(2:2:end)));
end
