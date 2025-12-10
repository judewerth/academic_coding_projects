function am = alfa_m(V)

am = .1 * (V+40) ./ (1 - exp(-(V+40)/10));

inds = find(isnan(am));
am(inds) = 1;