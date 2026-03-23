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
