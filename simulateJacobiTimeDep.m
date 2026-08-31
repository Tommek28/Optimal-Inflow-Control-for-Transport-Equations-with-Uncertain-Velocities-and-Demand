function [t, X] = simulateJacobiTimeDep(d0, kappa, theta_fun, sigma, dt, T, N)
% Simuliert N Jacobi-Prozesse auf [1,3] mit zeitabhängigem theta(t)
%
% INPUT:
% d0         - Startwert (Skalar oder 1xN Vektor)
% kappa      - Mean-Reversion Geschwindigkeit
% theta_fun  - Funktionshandle, z.B. @(t) 2 + 0.5*sin(t)
% sigma      - Diffusionsparameter
% dt         - Zeitschritt
% T          - Endzeit
% N          - Anzahl Pfade
%
% OUTPUT:
% t          - Zeitgitter (1xM)
% X          - Matrix (N x M), jeder Zeile = ein Pfad

    a = 0;
    b = 4;

    t = 0:dt:T;
    M = length(t);

    % Initialisierung
    X = zeros(N, M);
    
    if isscalar(d0)
        X(:,1) = d0;
    else
        X(:,1) = d0(:);
    end

    for i = 1:M-1
        dW = sqrt(dt) * randn(N, 1);   % N Zufallszahlen

        theta_t = theta_fun(t(i));     % Skalar
        
        drift = kappa * (theta_t - X(:,i));
        diffusion = sigma * sqrt(max((X(:,i)-a).*(b-X(:,i)), 0));

        X(:,i+1) = X(:,i) + drift * dt + diffusion .* dW;

        % Begrenzung auf [1,3]
        X(:,i+1) = min(max(X(:,i+1), a), b);
    end
end