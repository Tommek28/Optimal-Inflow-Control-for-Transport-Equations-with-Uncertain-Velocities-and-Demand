

function [outflow,rho] = transportUpwind(lambda, u, dt, T)
% Simuliert rho_t + lambda rho_x = 0 mit CFL=1 (Upwind)
%
% INPUT:
% lambda  - (N x 1) Geschwindigkeiten
% u       - Funktionshandle @(t) oder Vektor (1 x Nt)
% dt      - Zeitschritt
% T       - Endzeit
%
% OUTPUT:
% outflow - (N x Nt) Ausfluss bei x=1

    t = 0:dt:T;
    Nt = length(t);
    N = length(lambda);

    outflow = zeros(N, Nt);

    for k = 1:N   % über alle Pfade
        lam = lambda(k);

        dx = lam * dt;
        Nx = ceil(1 / dx);      % Anzahl Raumzellen
        dx = 1 / Nx;            % leicht angepasst

        rho = zeros(1, Nx);     % aktuelle Lösung

        for n = 1:Nt-1
    % Randwert
            if isa(u, 'function_handle')
                u_val = u(t(n));
            else
                u_val = u(n);
            end

            rho(1) = u_val / lam;

            % Schutz gegen unendliche / NaN
            if ~isfinite(rho(1))
                rho(1) = 0;
            end

            % Shift-Upwind
            rho(2:end) = rho(1:end-1);

            % Ausfluss
            outflow(k,n) = rho(end) * lam;
        end

        % letzter Zeitpunkt
        outflow(k, Nt) = rho(end)*lam;
    end
end
