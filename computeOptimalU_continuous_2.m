

function u_vals = computeOptimalU_continuous_2(t_grid, d0, kappa, theta_fun, alpha, beta, dt_int,proxy,velocityShift)
    N = length(t_grid) - 1;
  %  u_vals = zeros(1, N);
    m_fun = @(t) jacobiMeanClosed(t, d0, kappa, theta_fun);
    m_fun_vals = m_fun(t_grid(1):dt_int:t_grid(end));
    dl = linspace(velocityShift(1),velocityShift(2),100);
%    deltaL = dl(2)-dl(1);
    varphi = @(l) 1./(velocityShift(2)-velocityShift(1))*(l>=velocityShift(1)).*(l<=velocityShift(2));
    
    % Berechnung stetige optimale Kontrolle
    s = t_grid(1):dt_int:t_grid(end);
    u_vals = zeros(1,length(s));
    q = zeros(length(s),length(dl));
    for k=1:length(s)
        if max(s(k)+1./dl)>t_grid(end)
            k=k;
        end
        q(k,:) = (varphi(dl).*((s(k)+1./dl)>=1/velocityShift(1)).*((s(k)+1./dl)<=t_grid(end)))';
        if sum(q(k,:))==0
            k=k;
        else
            u_vals(k) = sum(m_fun(s(k)+1./dl).*q(k,:))/sum(q(k,:));
        end
    
    end
end