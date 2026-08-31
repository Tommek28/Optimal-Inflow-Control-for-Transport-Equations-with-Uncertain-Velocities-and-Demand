function prob = compute_probabilities(t_hat, t_star, s_grid, alpha, beta, c, N)
% Computes probabilities P_i(s_k) for all i,k
%
% INPUT:
% t_hat  : [0,...,1] injection grid (size I+1)
% t_star : velocity grid (size J+1)
% s_grid : evaluation times (size K+1)
% alpha, beta, c : parameters of shifted Beta
% N      : number of Monte Carlo samples
%
% OUTPUT:
% prob   : (K+1) x I matrix

I = length(t_hat) - 1;
J = length(t_star) - 1;
K = length(s_grid) - 1;

prob = zeros(K+1,I);

% --- Precompute lambda samples ---
lambda = c + betarnd(alpha, beta, N, J);

% --- Loop over s ---
for k = 1:(K+1)
    s = s_grid(k);

    % Precompute all weights for this s
    w_all = zeros(I+1, J);

    for i = 1:(I+1)
        t0 = t_hat(i);

        % compute weights w_j = overlap length
        left  = max(t0, t_star(1:J));
        right = min(s,  t_star(2:J+1));
        w = max(0, right - left);

        w_all(i,:) = w;
    end

    % --- Compute sums for all i in one go ---
    % S_i = sum_j w_{i,j} lambda_j
    S = lambda * w_all';   % size: N x (I+1)

    % --- Compute probabilities ---
    for i = 1:I
        S_i  = S(:,i);
        S_ip = S(:,i+1);

        prob(k,i) = mean( (S_i >= 1) & (S_ip < 1) );
    end
end

prob=prob';
end