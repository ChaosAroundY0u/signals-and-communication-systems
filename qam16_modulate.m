function qam16_out = qam16_modulate(x)
    qam16_out = 1./sqrt(10) .* ( (1-2.*x(1:4:end)) .* (2-(1-2.*x(3:4:end))) + 1j .* (1-2.*x(2:4:end)) .* (2-(1-2*x(4:4:end))) );
end
