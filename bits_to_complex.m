function out = bits_to_complex(symbs_count, nfft, num_subcarrier, modulation_technique)
    scale_var =                                     modulation_technique;
    x_count =                                       nfft * symbs_count * scale_var; % число бит
    x =                                             randi([0, 1], 1, x_count);
    out =                                           zeros(1, symbs_count .* nfft);
    switch scale_var
        case 1
            modulate =                              @bpsk_modulate;
        case 2
            modulate =                              @qpsk_modulate;
        case 4
            modulate =                              @qam16_modulate;
        case 6
            modulate =                              @qam64_modulate;
    end
    for iterable = 1:symbs_count
        x_use =                                     x(iterable*scale_var*nfft - scale_var*nfft + scale_var*num_subcarrier + 1:iterable*scale_var*nfft - scale_var * num_subcarrier);
        %symbol =                                    x(iterable*scale_var*nfft-scale_var*nfft+1:iterable*scale_var*nfft);
        modulated =                                 modulate(x_use);
        out(iterable*nfft-nfft+1:iterable*nfft) =   [zeros(1, num_subcarrier), modulated, zeros(1, num_subcarrier)];
    end
end
