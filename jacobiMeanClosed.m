


function m_t = jacobiMeanClosed(t_query, d0, kappa, theta_fun)
% Geschlossene Darstellung mit numerischer Quadratur
%
% m(t) = e^{-kappa t} d0 + kappa * int_0^t e^{-kappa(t-s)} theta(s) ds

    m_t = zeros(size(t_query));

    for j = 1:length(t_query)
        t = t_query(j);

        % Parameter
        Ns = 1000;                 % Anzahl Stützstellen (fein genug wählen!)
        s = linspace(0, t, Ns);    % Gitter
        
        ds = s(2) - s(1);          % Schrittweite
        
        % Integrand auswerten
        vals = exp(-kappa*(t - s)) .* theta_fun(s);
        
        % Rechteckregel (links)
        I = sum(vals(1:end-1) + vals(2:end))/2 * ds;

        % integrand = @(s) exp(-kappa*(t - s)) .* theta_fun(s);
        % 
        % I = integral(integrand, 0, t);   % adaptive Quadratur

        m_t(j) = exp(-kappa*t)*d0 + kappa * I;
    end
end