function bits = hard_demapping(symbols)
     symbols = reshape(symbols.', 1, []);
     constellation = [1+1i, 1-1i, -1+1i, -1-1i] / sqrt(2);
     bits_map = [0 0; 0 1; 1 0; 1 1];
     bits = [];
     for sym = symbols
        % Находим ближайшую точку созвездия
        [~, idx] = min(abs(sym - constellation));
        bits = [bits, bits_map(idx, :)];
     end
 end
