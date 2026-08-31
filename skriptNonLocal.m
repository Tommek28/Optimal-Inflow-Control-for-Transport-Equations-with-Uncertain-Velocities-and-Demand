% Ein Skript für das nonlocal-Modell von Timo

xrange = [0,1];
uFct = @(t) 0.2 + 0.05*sin(t);

dx = 0.05;
dt = dx/2;
T = 10;
u = uFct(0:dt:T);
v = @(rho) 1-rho.^2;
eta = 0.2;
W = @(x) 3/(2*eta)*(eta^2-x.^2);

%Zufall einbauen, Störung um max Amplitude 0.1 nach oben und unten
MC = 100;
epsilon = 0.4*rand(1,MC) - 0.2;
outflow = zeros(T/dt,MC);
tic()
for i=1:MC
    v = @(rho) max(0,1-rho.^2 + epsilon(i));
    [outflow(:,i),~] = stochasticNonlocalModel(xrange,u,dx,dt,T,v,W,eta);
end
tid=toc()


%[outflow,rho] = stochasticNonlocalModel(xrange,u,dx,dt,T,v,W,eta);

figure
plot(dt:dt:T, mean(outflow,2), LineWidth=2)
hold on
plot(dt:dt:T, median(outflow,2), LineWidth=1.5, LineStyle="--")
plot(dt:dt:T, quantile(outflow,0.2,2), LineWidth=1.2, LineStyle=":")
plot(dt:dt:T, quantile(outflow,0.8,2), LineWidth=1.2, LineStyle="-.")
xlabel('time')
ylabel('outflow')
legend('Mean Outflow','Median Outflow', '20% quantile', '80% quantile','Location','best')
title('Outflow of nonlocal model')

% figure
% for i=1:10:(T/dt)
%     plot((xrange(1)+dx):dx:(xrange(2)-dx),rho(i,1:end-1),'LineWidth',2)
%     ylim([0,1])
%     pause(0.0001)
% end


%% Demand und optimale Kontrolle für das deterministische Problem - Spielwiese
% Optimierungszeitraum ist [3,T], die ersten drei Zeiteinheiten werden zum
% Einschwingen genutzt

D = @(t) 0.1 + 0.05*cos(t);
DD = D(3:dt:T)';


d0 = 0.16;
kappa=4;
theta_fun = @(t) 0.1 + 0.05*cos(t);
D = jacobiMeanClosed(0:dt:T, d0, kappa, theta_fun);

DD = D(3/dt+1:end)';



v = @(rho) 1-rho.^2;
objFct = @(u) dt*sum((DD-stochasticNonlocalModel_control(xrange,u,dx,dt,T,v,W,eta)).^2);

tic()
[optInflow, residuum] = fmincon(objFct,0.2*ones(1,length(0:dt:T)));

tidOpt = toc()



figure()
plot(0:dt:T,optInflow, 'LineWidth',1.2,'LineStyle',':')
hold on
plot(dt:dt:T,stochasticNonlocalModel(xrange,optInflow,dx,dt,T,v,W,eta),'LineWidth',2)
plot(0:dt:T, D(0:dt:T), 'LineWidth',1.5,'LineStyle','--')
hold off
legend('Inflow','outflow','demand')
xlabel('time')
ylabel('inflow')
title('Comparison of optimal outflow and demand - proxy')


%% Jetzt der echte Proxy ... und danach das stochastische Problem

d0 = 0.16;
kappa=4;
theta_fun = @(t) 0.1 + 0.05*cos(t);
D = jacobiMeanClosed(0:dt:T, d0, kappa, theta_fun);

DD = D(3/dt+1:end)';

%D = @(t) 0.1 + 0.05*cos(t);
%DD = D(3:dt:T)';


v = @(rho) ...
    (1 - rho.^2).*(1 - rho.^2 >= 0.2) + ...
    ((1 - rho.^2) + 1.25*(0.2 - (1 - rho.^2)).^2).*(1 - rho.^2 < 0.2);
disp('Achtung diese Geschwindigkeit gilt nur für U([-0.2,0.2])!! das wurde unten manchmal geändert')
objFct = @(u) dt*sum((DD-stochasticNonlocalModel_control(xrange,u,dx,dt,T,v,W,eta)).^2);

tic()
[optInflowProxy, residuumProxy] = fmincon(objFct,0.15*ones(1,length(0:dt:T)));

tidOptProxy = toc()



figure()
plot(0:dt:T,optInflowProxy, 'LineWidth',1.2,'LineStyle',':')
hold on
plot(dt:dt:T,stochasticNonlocalModel(xrange,optInflowProxy,dx,dt,T,v,W,eta),'LineWidth',2)
plot(0:dt:T, D, 'LineWidth',1.5,'LineStyle','--')
hold off
legend('Inflow','outflow','demand')
xlabel('time')
ylabel('inflow')
title('Comparison of optimal outflow and demand - proxy')


% Optimale Kontrolle für das stochastische Problem
% wir müssen hierfür 
% - zuerst eine MC-Simulation für die ersten beiden Ausflussmomente
% - Danach optimierung wie vorher

MC = 200;
v_vec = cell(MC,1);
epsilon = 0.4*rand(1,MC) - 0.2;
for i=1:MC
    v_vec{i} = @(rho) max(0,1-rho.^2 + epsilon(i));
end

%[output] = stochasticNonlocalModel_control_stochastic(xrange,u,dx,dt,T,v,W,eta,DD);

objFctStoch = @(u)stochasticNonlocalModel_control_stochastic(xrange,u,dx,dt,T,v_vec,W,eta,DD) ;


tic()
[optInflow_Stoch, residuum_Stoch] = fmincon(objFctStoch,0.15*ones(1,length(0:dt:T)));

tidOptStoch = toc()


forward_Out_Stoch = stochasticNonlocalModel_forward(xrange,optInflow_Stoch,dx,dt,T,v_vec,W,eta);
forward_Out_Proxy = stochasticNonlocalModel_forward(xrange,optInflowProxy,dx,dt,T,v_vec,W,eta);



figure()
plot(0:dt:T,optInflow_Stoch, 'LineWidth',1,'LineStyle','-','Color','r')
hold on
plot(0:dt:T,optInflowProxy, 'LineWidth',1,'LineStyle','-','Color','b')
%plot(dt:dt:T,stochasticNonlocalModel(xrange,optInflow_Stoch,dx,dt,T,v,W,eta),'LineWidth',2,'Color','r')
plot(dt:dt:T,mean(forward_Out_Stoch,2),'LineWidth',2,'Color','r','LineStyle',':')
%plot(dt:dt:T,stochasticNonlocalModel(xrange,optInflowProxy,dx,dt,T,v,W,eta),'LineWidth',2,'Color','b')
plot(dt:dt:T,mean(forward_Out_Proxy,2),'LineWidth',2,'Color','b','LineStyle',':')
plot(0:dt:T, D, 'LineWidth',1.5,'LineStyle','--','Color','g')
hold off
legend('Inflow Stochastic','Inflow Proxy','Outflow Stochastic','Outflow Proxy','Demand','Location','best')
xlabel('time')
ylabel('inflow')
title('Comparison of optimal outflow and demand - Stochastic vs. Proxy')

% Fehler in der Zielfunktion
disp('Error Stochastic')
err_Stoch = dt*sum(mean((D(3/dt+1:end)' - forward_Out_Stoch(3/dt:end,:)).^2,2))
disp('Error Proxy')
err_Proxy = dt*sum(mean((D(3/dt+1:end)'- forward_Out_Proxy(3/dt:end,:)).^2,2))

% Test

% disp('Res stochastic')
% stochasticNonlocalModel_control_stochastic(xrange,optInflowProxy,dx,dt,T,v_vec,W,eta,DD)
% disp('Res det Proxy')
% stochasticNonlocalModel_control_stochastic(xrange,optInflow_Stoch,dx,dt,T,v_vec,W,eta,DD)
% 
% figure
% plot(dt:dt:T,stochasticNonlocalModel(xrange,optInflowProxy,dx,dt,T,v,W,eta),'LineWidth',2)
% 
% hold on
% plot(dt:dt:T,stochasticNonlocalModel(xrange,optInflow_Stoch,dx,dt,T,v,W,eta),'LineWidth',2)
% plot(0:dt:T,D(0:dt:T))
% hold off
% legend('Det','Stochastisch','Demand')


%% Eine Varianzuntersuchung

k=-5:-1;

meanErrorStoch = zeros(length(k),1);
meanErrorProxy = zeros(length(k),1);

v_helpProxy = @(rho,tau) ...
    ( (1 - rho.^2 >= tau) .* (1 - rho.^2) ) + ...
    ( (1 - rho.^2 >= 0 & 1 - rho.^2 < tau) .* ...
      ((1 - rho.^2) + (tau - (1 - rho.^2)).^2 ./ (4*tau)) ); 
idx = 0;
for qqq = k
    qqq
    idx=idx+1;
   % D = @(t) 0.1 + 0.05*cos(t);
    %DD = D(3:dt:T)';
   % v = @(rho) ...
   %    (1 - rho.^2).*(1 - rho.^2 >= 0.2) + ...
   %    ((1 - rho.^2) + 1.25*(0.2 - (1 - rho.^2)).^2).*(1 - rho.^2 < 0.2);
   
    v = @(rho) v_helpProxy(rho,2^(qqq));
   disp('Achtung diese Geschwindigkeit gilt nur für U([-0.2,0.2])!! das wurde unten manchmal geändert')
    objFct = @(u) dt*sum((DD-stochasticNonlocalModel_control(xrange,u,dx,dt,T,v,W,eta)).^2);
    
    tic()
    [optInflowProxy, residuumProxy] = fmincon(objFct,0.15*ones(1,length(0:dt:T)));
    
    tidOptProxy = toc()
    
    
    % Optimale Kontrolle für das stochastische Problem
    % wir müssen hierfür 
    % - zuerst eine MC-Simulation für die ersten beiden Ausflussmomente
    % - Danach optimierung wie vorher
    
    MC = 200;
    v_vec = cell(MC,1);
    %%epsilon = 0.4*rand(1,MC) - 0.2;
    epsilon = 2*(2^(qqq))*rand(1,MC)-2^(qqq);
    for i=1:MC
        v_vec{i} = @(rho) max(0,1-rho.^2 + epsilon(i));
    end
    
    %[output] = stochasticNonlocalModel_control_stochastic(xrange,u,dx,dt,T,v,W,eta,DD);
    
    objFctStoch = @(u)stochasticNonlocalModel_control_stochastic(xrange,u,dx,dt,T,v_vec,W,eta,DD) ;
    
    
    tic()
    [optInflow_Stoch, residuum_Stoch] = fmincon(objFctStoch,0.15*ones(1,length(0:dt:T)));
    
    tidOptStoch = toc()
    
    
    forward_Out_Stoch = stochasticNonlocalModel_forward(xrange,optInflow_Stoch,dx,dt,T,v_vec,W,eta);
    forward_Out_Proxy = stochasticNonlocalModel_forward(xrange,optInflowProxy,dx,dt,T,v_vec,W,eta);

    meanErrorStoch(idx) = dt*sum(mean((D(3/dt+1:end)' - forward_Out_Stoch(3/dt:end,:)).^2,2));

    meanErrorProxy(idx) = dt*sum(mean((D(3/dt+1:end)' - forward_Out_Proxy(3/dt:end,:)).^2,2));

end

figure
plot(k, log2(meanErrorStoch), 'LineWidth',2)
hold on
plot(k, log2(meanErrorProxy), 'LineWidth',2)
hold off
xlabel('k')
ylabel('log2(error)')
legend('Stochastic','Proxy')
title('Error evolution for smaller supports of \varepsilon')

figure
plot(k, log2(abs(meanErrorStoch-meanErrorProxy)), 'LineWidth',2)
xlabel('k')
ylabel('log2(error)')
legend('Excess cost')
title('Error evolution for smaller supports of \varepsilon in the nonlocal model')


figure

y = log2(abs(meanErrorStoch - meanErrorProxy));

y=y(1:end-1)
k1 =k(1:end-1)

% Daten plotten
plot(k1, y, 'o-', 'LineWidth', 2, 'MarkerSize', 7)
hold on

% Linearer Best Fit
p = polyfit(k1, y, 1);
yFit = polyval(p, k1);

plot(k1, yFit, '--', 'LineWidth', 2)

% Steigung
slope = p(1);

xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\log_2(\mathrm{error})$', 'Interpreter', 'latex')

legend( ...
    'Excess cost', ...
    sprintf('Linear fit, slope = %.3f', slope), ...
    'Location', 'best' ...
)

title('Error evolution for smaller supports of $\varepsilon$ in the nonlocal model', ...
    'Interpreter', 'latex')

grid on
box on