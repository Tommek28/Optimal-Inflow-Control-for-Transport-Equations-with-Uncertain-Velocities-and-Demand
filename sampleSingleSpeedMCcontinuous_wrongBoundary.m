function [L1_error_mean,L2_error_mean,L3_error_mean,mean_out,mean_X,lambda,X,u]  = sampleSingleSpeedMCcontinuous_wrongBoundary(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,ii,velocityShift,proxy,randeinzug)

% Jacobi-Nachfrage 
% Beta Verteilung auf der Kante mit einer fixen Geschwindigkeit
% Rückgabe: durchschnittlicher L^1, L^2, L^3 Fehler, durchschnittlicher outflow,
% durchschnittliche Demand

t_grid = 0:dt:T;

% Sample Jacobi Nachfrage
[t, X] = simulateJacobiTimeDep(d0, kappa, theta_fun, sigma, dt, T,MCruns);

% Geschwindigkeit
lambda = velocityShift(1)+(velocityShift(2)-velocityShift(1))*betarnd(alpha,beta,MCruns,1);
%lambda = 0.5+ones(MCruns,1);

t_gridU = 0:dt:T;

% falls proxy dann wird hier die Geschwindigkeitsverteilung nahe ein dirac
% delta!
if proxy
    u = computeOptimalU_Riemann_2_wrongBoundary(t_gridU, d0, kappa, theta_fun, alpha, beta, 5*10^-4,0,[2-1000*eps,2+1000*eps]);
else
    u = computeOptimalU_Riemann_2_wrongBoundary(t_gridU, d0, kappa, theta_fun, alpha, beta, 5*10^-4,0,velocityShift);
end

% u = [u1, u2, ..., u10]   % stückweise konstant
u_t = zeros(size(t_grid));

N = length(u);  % Anzahl Intervalle

for i = 1:N-1
    idx = (t >= t_gridU(i)) & (t < t_gridU(i+1));
    u_t(idx) = u(i);
end
% Optional: letzter Punkt T
u_t(t_grid >= t_gridU(end)) = u(end);

% Ausflussberechnung mit Upwind
out = transportUpwind(lambda, u_t, dt, T);


% Plotting und Fehler
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
        plot(t, out_k, 'Color', colors(k,:), 'LineWidth', 1.5);
    
        % Plot Jacobi-Nachfrage (dotted)
        plot(t, X_k, '--', 'Color', colors(k,:), 'LineWidth', 1);

        
    end


    xlabel('t');
    ylabel('Outflow / Demand');
    title('Outflow vs. Jacobi-Demand (Realizations)');
    legend_entries = cell(1,min(MCruns,1)*2);
    for k = 1:2:(2*min(MCruns,1))
        legend_entries{k} = ['MC Outflow ' num2str(k)];
    end
    for k = 2:2:(2*min(MCruns,1))
        legend_entries{k} = ['MC Demand ' num2str(k)];
    end
    legend(legend_entries);
    hold off;

    figure(2*ii);
    hold on;
    mean_out = mean(out);
    mean_X = mean(X);
    % Mittelwert Ausfluss (solid, blau)
    plot(t, mean_out, '-b', 'LineWidth', 2);
    
    % Mittelwert Jacobi-Nachfrage (dotted, rot)
    plot(t, mean_X, '--r', 'LineWidth', 2);
    
    xlabel('t');
    ylabel('Outflow / Demand');
    title('Outflow vs. Jacobi-Demand (Realizations)');
    legend('Mean Outflow','Mean Demand');
    grid on;
    hold off;

end

% Anzahl Zeitschritte
N = length(t);

% Fehler pro Pfad ab t=1!!!

tminIdx = randeinzug/dt+1;
tmaxIdx = (T-randeinzug)/dt+1;
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



% Mittelwert über alle MC-Pfade
mean_out = mean(out, 1);  % 1 x N
mean_X   = mean(X, 1);    % 1 x N

