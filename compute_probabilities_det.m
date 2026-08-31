function prob = compute_probabilities_det(t_hat, t_star, s_grid, lambda)

I = length(t_hat) - 1;
J = length(t_star) - 1;
K = length(s_grid) - 1;

prob = zeros(K+1, I);

for k = 1:(K+1)
    s = s_grid(k);

    % Gewichte für alle i gleichzeitig
    w_all = zeros(I+1, J);

    for i = 1:(I+1)
        t0 = t_hat(i);

        left  = max(t0, t_star(1:J));
        right = min(s,  t_star(2:J+1));

        w_all(i,:) = max(0, right - left);
    end

    % deterministische Summen
    S = w_all * lambda(:);   % (I+1) x 1

    % Indikatorbedingung
    for i = 1:I
        prob(k,i) = (S(i) >= 1) && (S(i+1) < 1);
    end
end
prob=prob';
end