

function u_vals = computeOptimalU_continuous(t_grid, d0, kappa, theta_fun, alpha, beta, dt_int,proxy,velocityShift)
    N = length(t_grid) - 1;
    u_vals = zeros(1, N);
    m_fun = @(t) jacobiMeanClosed(t, d0, kappa, theta_fun);
    m_fun_vals = m_fun(t_grid(1):dt_int:t_grid(end));
    q=diff(t_grid);
        for i = 1:N
            
            ti = t_grid(i);
            ti1 = t_grid(i+1);
    
            % s-Gitter von ti bis T
            s = (ti+dt_int):dt_int:t_grid(end);
            if proxy
                u_vals(i) = m_fun(ti + (alpha+beta)/(alpha*velocityShift(2)+beta*velocityShift(1)));
                
            else
                %w_vals = 1./(s-ti).^2.*betapdf((1./(s-ti) -velocityShift(1))/(velocityShift(2)-velocityShift(1)), alpha, beta);
                w_vals = 1./(s-ti).^2.*unifpdf((1./(s-ti)), velocityShift(1), velocityShift(2));
                m_vals = m_fun_vals(end-length(w_vals)+1:end);%m_fun(s);
        
                ds = diff([s, s(end)+dt_int]);  % letzter Schritt dt_int
                num = sum(m_vals .* w_vals .* ds);
                den = sum(w_vals .* ds);
                
                u_vals(i) = num / den;
                %u_vals(i) = mean(m_fun_vals(  ((i+1)+round(1/(3*dt_int))):min(length(m_fun_vals),((i+1)+round(1/dt_int))) ));
                %u_vals(i) = mean(m_fun_vals(  )   )
            end
        end
end