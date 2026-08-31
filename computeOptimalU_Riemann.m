

function u_vals = computeOptimalU_Riemann(t_grid, d0, kappa, theta_fun, alpha, beta, dt_int,proxy,velocityShift)
    N = length(t_grid) - 1;
    u_vals = zeros(1, N);
    m_fun = @(t) jacobiMeanClosed(t, d0, kappa, theta_fun);
    m_fun_vals = m_fun(t_grid(1):dt_int:t_grid(end));
    q=diff(t_grid);
        for i = 1:N
            ti = t_grid(i);
            ti1 = t_grid(i+1);
    
            % s-Gitter von ti bis T
            s = ti:dt_int:t_grid(end);

            w_vals = weightFunction_Riemann(s, ti, ti1, alpha, beta,proxy,velocityShift);
           % idx = max(1, round(1 + (i-1)*q(1)/dt_int));

            m_vals = m_fun_vals(end-length(w_vals)+1:end);%m_fun(s);
    
            ds = diff([s, s(end)+dt_int]);  % letzter Schritt dt_int
            num = sum(m_vals .* w_vals .* ds);
            den = sum(w_vals .* ds);
            
            
            u_vals(i) = num / den;
        end
end