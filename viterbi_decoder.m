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
