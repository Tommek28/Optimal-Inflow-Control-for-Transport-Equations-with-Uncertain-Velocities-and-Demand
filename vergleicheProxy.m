function [L1_error_mean,L2_error_mean,L3_error_mean,mean_out,mean_X,u] = vergleicheProxy(lambda,X,dt,T,dtU,d0,kappa,theta_fun,alpha,beta,MCruns,ii,plotBool,velocityShift,randeinzug)
t = 0:dt:T;
t_gridU = 0:dtU:T;
% Berechne optimalen Inflow

u = computeOptimalU_Riemann_2(t_gridU, d0, kappa, theta_fun, alpha, beta, 5*10^-4,1,[2-1000*eps,2+1000*eps]);
u_t = zeros(size(t));

N = length(u);  % Anzahl Intervalle

for i = 1:N
    idx = (t >= t_gridU(i)) & (t < t_gridU(i+1));
    u_t(idx) = u(i);
end
% Optional: letzter Punkt T
u_t(t >= t_gridU(end)) = u(end);
out = transportUpwind(lambda, u_t, dt, T);


tminIdx = randeinzug/dt+1;
tmaxIdx = (T-randeinzug)/dt +1 ;
err2 = zeros(MCruns,1);
err1 = zeros(MCruns,1);
err3 = zeros(MCruns,1);
Xerr = X(:,tminIdx:tmaxIdx);
outErr = out(:,tminIdx:tmaxIdx);
for k = 1:MCruns
    % Trapezregel über das Zeitgitter t
    err1(k) = trapz(t(tminIdx:tmaxIdx), abs(Xerr(k,:) - outErr(k,:)));  % L1-Fehler
    err2(k) = trapz(t(tminIdx:tmaxIdx), (Xerr(k,:) - outErr(k,:)).^2);  % L2-Fehler^2
    err3(k) = trapz(t(tminIdx:tmaxIdx), (Xerr(k,:) - outErr(k,:)).^3);  % L3-Fehler^3
end

% Mittelwert über alle MC-Pfade
L1_error_mean = mean(err1);
L2_error_mean = mean(err2);
L3_error_mean = mean(err3);

mean_out = mean(out);
mean_X = mean(X);

colors = lines(MCruns);  % MATLAB-Farbpalette
if plotBool
    figure(2*ii-1);
    hold on;
    
    for k = 1:min(MCruns,1)
        % Ausfluss Pfad k
        out_k = out(k,:);
    
        % Jacobi-Pfad k
        X_k = X(k,:);
    
        % Plot Ausfluss (solid)
        plot(t, out_k, '-.' , 'Color', colors(k+10,:), 'LineWidth', 1.5);
    
        % Plot Jacobi-Nachfrage (dotted)
        %plot(t, X_k, '--', 'Color', colors(k,:), 'LineWidth', 1);

        
    end


    % xlabel('t');
    % ylabel('Outflow / Demand');
    % title('Outflow vs. Jacobi-Demand (Realizations)');
    % legend_entries = cell(1,min(MCruns,3)*2);
    % for k = 7:1:9
    %     legend_entries{k} = ['Proxy Outflow ' num2str(k)];
    % end
    % 
    % legend(legend_entries);
    hold off;

    figure(2*ii);
    hold on;
    mean_out = mean(out);
    mean_X = mean(X);
    % Mittelwert Ausfluss (solid, blau)
    plot(t, mean_out, ':g', 'LineWidth', 2);
    
   
    
   % xlabel('t');
    %ylabel('Outflow / Demand');
    %title('Outflow vs. Jacobi-Demand (Realizations)');
    %legend('Mean Outflow','Mean Demand');
    grid on;
    hold off;

end



end