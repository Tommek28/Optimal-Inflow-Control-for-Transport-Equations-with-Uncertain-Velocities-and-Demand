    % nur für eine gleichverteilung nutzen

function u_vals = computeOptimalU_Riemann_2_wrongBoundary(t_grid, d0, kappa, theta_fun, alpha, beta, dt_int,proxy,velocityShift)
N = length(t_grid) - 1;
u_vals = zeros(1, N);
m_fun = @(t) jacobiMeanClosed(t, d0, kappa, theta_fun);
%    m_fun_vals = m_fun(t_grid(1):dt_int:t_grid(end));
q=diff(t_grid);
dl = linspace(velocityShift(1),velocityShift(2),100);
%    deltaL = dl(2)-dl(1);
varphi = @(l) 1./(velocityShift(2)-velocityShift(1))*(l>=velocityShift(1)).*(l<=velocityShift(2));

% Berechnung stetige optimale Kontrolle
s = t_grid(1):dt_int:t_grid(end);
uStar = zeros(1,length(s));
q = zeros(length(s),length(dl));
for k=1:length(s)
    q(k,:) = (varphi(dl))';             % Hier liegt die einzige Änderung im Vergleich zu vorher
    if sum(q(k,:))==0
        k=k;
    else
        uStar(k) = sum(m_fun(s(k)+1./dl).*q(k,:))/sum(q(k,:));
    end

end


    for i = 1:N
        ti = t_grid(i);
        ti1 = t_grid(i+1);

        % s-Gitter von ti bis T
        s = ti:dt_int:t_grid(end);

       %  w_vals = weightFunction_Riemann(s, ti, ti1, alpha, beta,proxy,velocityShift);
       % % idx = max(1, round(1 + (i-1)*q(1)/dt_int));
       % 
       %  m_vals = m_fun_vals(end-length(w_vals)+1:end);%m_fun(s);
       % 
       %  ds = diff([s, s(end)+dt_int]);  % letzter Schritt dt_int
       %  num = sum(m_vals .* w_vals .* ds);
       %  den = sum(w_vals .* ds);
        
        
     %   u_vals(i) = num / den;
     idx1 = round(t_grid(i)/dt_int) + 1;
     idx2 = round(t_grid(i+1)/dt_int) + 1;
     u_vals(i) = sum(sum(q(idx1:idx2,:),2)'.*uStar(idx1:idx2))/sum(sum(q(idx1:idx2,:),2));

    end
end