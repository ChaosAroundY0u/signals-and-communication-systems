function v = dec2binvec(d, n)
    v = mod(floor(d ./ 2.^(0:n-1)), 2);
end
