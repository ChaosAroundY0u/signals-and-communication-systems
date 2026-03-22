clear;
close all;
a = [0 1 0 1 1 1 1 0 0 0 1 ];
%x_count = 4095;
%x = randi([0, 1], 1, x_count);
%коэффициенты полиномов идут слева направо, т.е.:
%g1 = [1 0]; % 0*x + 1 или  вернее 1 + 0*x
%g2 = [1 1]; % 1*x + 1 -/-/-/-/-/- 1 + 1*x
g1 = [1 1 0 1 1 0 1]; % 133 8
g2 = [1 0 0 1 1 1 1]; % 171 8
%под знаком сложения подразумевается операция xor (сложение по модулю 2)

encoded = conv_encoder(a, g1, g2);
decoded = viterbi_decoder(encoded, g1, g2);

function y = conv_encoder(x, g1, g2)
    % rate 1/2
    n1 =        size(g1);
    n2 =        size(g2);
    y = [];
    if ~isequal(n1, n2)
        return;
    end
    reg =       zeros(n1);
    in_size =   size(x);
    nx =        max(in_size);
    r =         max(n1);
    out_size =  2*nx + (r-1)*2;
    if (in_size(1) == 1)
        y = zeros(1, out_size);
    elseif (in_size(2) == 1)
        y = zeros(out_size, 1);
    else
        return;
    end

    for i = 1:nx
        reg(1) = x(i);
        y(2*i-1) = mod(sum(reg.*g1), 2);
        y(2*i)   = mod(sum(reg.*g2), 2);
        for k = 2:r
            reg(r-k+2) = reg(r-k+1);
        end
    end

    for i = nx+1:nx+r-1
        reg(1) = 0;
        y(2*i-1) = mod(sum(reg.*g1), 2);
        y(2*i)   = mod(sum(reg.*g2), 2);
        for k = 2:r
            reg(r-k+2) = reg(r-k+1);
        end
    end
end

function x_hat = viterbi_decoder(y, g1, g2)
    r = length(g1);
    num_states = 2^(r-1);
    nx = (length(y) / 2) - (r - 1);
    %метрики путей
    PM = inf(num_states, 1);
    PM(1) = 0;

    paths = zeros(num_states, nx + r - 1);

    for i = 1:nx + r - 1
        new_PM = inf(num_states, 1);
        new_paths = zeros(num_states, nx + r - 1);
        current_bits = y(2*i-1 : 2*i);

        for state = 0:num_states-1
            if PM(state + 1) == inf, continue;
            end

            for bit = 0:1
                reg = [bit, dec2binvec(state, r-1)];

                out1 = mod(sum(reg .* g1), 2);
                out2 = mod(sum(reg .* g2), 2);

                dist = sum(mod(current_bits + [out1, out2], 2)); % длина метрики

                next_state = bit * 2^(r-2) + floor(state / 2);

                if PM(state + 1) + dist < new_PM(next_state + 1)
                    new_PM(next_state + 1) = PM(state + 1) + dist;
                    new_paths(next_state + 1, :) = paths(state + 1, :);
                    new_paths(next_state + 1, i) = bit;
                end
            end
        end
        PM = new_PM;
        paths = new_paths;
    end

    best_path = paths(1, :);
    x_hat = best_path(1:nx);
end

function v = dec2binvec(d, n)
    v = mod(floor(d ./ 2.^(0:n-1)), 2);
end
