function an = alfa_n(V)

an = .01 * (V+55) ./ (1 - exp(-(V+55)/10));

inds = find(isnan(an));
an(inds) = 0.1;