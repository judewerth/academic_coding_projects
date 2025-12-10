function bh = beta_h(V)

bh = 1 ./ (1 + exp(-(V+35)/10));