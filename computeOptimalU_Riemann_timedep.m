

function u_vals = computeOptimalU_Riemann_timedep(t_grid, d0, kappa, theta_fun, alpha, beta, dt_int,proxy,lambdaMin,lambdaMax,tiHat,velocityShift)
    N = length(t_grid) - 1;     % das ist das u-Gitter
    u_vals = zeros(1, N);
    m_fun = @(t) jacobiMeanClosed(t, d0, kappa, theta_fun);
    m_fun_vals = m_fun(t_grid(1):dt_int:t_grid(end));
    q=diff(t_grid);
    
    s_grid = 0:dt_int:t_grid(end);
    NMC = 2*10^5;

    
    if proxy
        prob = compute_probabilities_det(t_grid, tiHat, s_grid, velocityShift+alpha/(alpha+beta)*ones(length(tiHat)-1,1));
    else
        %prob = plj_diffGrids(alpha,beta,tiHat,0:dt_int:t_grid(end),lambdaMin,lambdaMax,t_grid,velocityShift);
        prob = compute_probabilities(t_grid, tiHat, s_grid, alpha, beta, velocityShift, NMC);
    end
    
    Elambda = (velocityShift+alpha/(alpha+beta))*ones(1,length(t_grid(1):dt_int:t_grid(end)));
    
        for i = 1:N
            ti = t_grid(i);
            ti1 = t_grid(i+1);
    
            % s-Gitter von ti bis T
            s = ti:dt_int:t_grid(end);
            
            
            w_vals = prob(i,((i-1)*t_grid(end)/N/dt_int+1):end); %weightFunction_Riemann(s, ti, ti1, alpha, beta,proxy);
            m_vals = m_fun_vals((1+(i-1)*q(1)/dt_int):end);%m_fun(s);
            
            ds = diff([s, s(end)+dt_int]);  % letzter Schritt dt_int
            num = nansum(Elambda(((i-1)*t_grid(end)/N/dt_int+1):end).*m_vals .* w_vals .* ds);
            den = nansum(Elambda(((i-1)*t_grid(end)/N/dt_int+1):end).*w_vals .* ds);
            
            
            u_vals(i) = num / den;
            
        end
    
end