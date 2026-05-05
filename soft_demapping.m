function bits = soft_demapping(symbols, noise_power, modulation_technique)
     scale_var = modulation_technique;
     noise_power_linear = 10.^(noise_power./10);
     switch scale_var
         case 2
             constellation = [1 1; 1 -1; -1 1; -1 -1];
             %bits_map =     [0 0; 0  1;  1 0;  1  1];
             bits =          zeros(1, scale_var*length(symbols));
             euDists =       scale_var .^ 2; 
             euclid_dist =   zeros(1, euDists);
                  for k = 1:length(symbols)
                    for m = 1:euDists
                        euclid_dist(1, m) = ((real(symbols(k)) - constellation(m, 1)).^2 + (imag(symbols(k)) - constellation(m, 2)).^2) ./ (2.*noise_power_linear);
                    end
                    LLR_b0 = log((exp(-euclid_dist(1)) + exp(-euclid_dist(2)))./(exp(-euclid_dist(3)) + exp(-euclid_dist(4))));
                    LLR_b1 = log((exp(-euclid_dist(1)) + exp(-euclid_dist(3)))./(exp(-euclid_dist(2)) + exp(-euclid_dist(4))));

                    bits(1, k*scale_var-scale_var+1:k*scale_var) = [LLR_b0, LLR_b1];
                  end
         case 4
             constellation = [1  1; 1  3; 3  1; 3  3; 1  -1; 1  -3; 3  -1; 3  -3; ...
                             -1  3;-1  1;-3  1;-3  3;-1  -1;-1  -3;-3  -1;-3  -3];
             %bits_map =     [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1;
             %                1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1; 1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1];
             euDists =       scale_var .^ 2;
             euclid_dist =   zeros(1, euDists);
             bits = zeros(1, scale_var*length(symbols));
             for k = 1:length(symbols)
                for m = 1:euDists
                    euclid_dist(1, m) = ((real(symbols(k)) - constellation(m, 1)).^2 + (imag(symbols(k)) - constellation(m, 2)).^2) ./ (2.*noise_power_linear);
                end
                LLR_b0 = log(sum(exp(-euclid_dist(1:8)))./(sum(exp(-euclid_dist(9:end)))));
                LLR_b1 = log((sum(exp(-euclid_dist(1:4))) + sum(exp(-euclid_dist(9:12)))) ./ (sum(exp(-euclid_dist(5:8))) + sum(exp(-euclid_dist(13:end)))));
                LLR_b2 = log((sum(exp(-euclid_dist(1:2))) + sum(exp(-euclid_dist(5:6))) + sum(exp(-euclid_dist(9:10))) + sum(exp(-euclid_dist(13:14)))) ./ (sum(exp(-euclid_dist(3:4))) + sum(exp(-euclid_dist(7:8))) + sum(exp(-euclid_dist(11:12))) + sum(exp(-euclid_dist(15:16)))));
                LLR_b3 = log(sum(exp(-euclid_dist(1:2:end))) ./ sum(exp(-euclid_dist(2:2:end))));

                bits(1, k*scale_var-scale_var+1:k*scale_var) = [LLR_b0, LLR_b1, LLR_b2, LLR_b3];
             end
             
     end
