function qam64_out = qam64_modulate(x)
    qam64_out = 1/sqrt(42) .* ( (1-2*x(1:6:end)) .* ( 4 - (1 - 2*x(3:6:end)).*(2 - (1 - 2*x(5:6:end))) ) + 1j * (1 - 2*x(2:6:end)) .* (4 - (1 - 2*x(4:6:end)) .* (2 - (1 - 2*x(6:6:end)))) );
end
