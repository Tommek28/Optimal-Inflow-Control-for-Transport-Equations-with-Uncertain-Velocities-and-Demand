
%% Content

% 1. Analyis von (t_{i+1}-t_i)->0 
% Transistion from piecewise continuous to continuous, including proxy

% 2. Analysis of varying the Variance of the transport velocity
% including proxy

% 3. Bonudary effects



%% Decreasing the length of the piecewise constant sections (figure 2,3)

% 1. Analyse von (t_{i+1}-t_i)->0 in den settings 1,3,4, also der Fehler im
% Übergang von piecewise continuous and continuous

% 2. Die Analyse hinsichtlich der Varianz in der Geschwindigkeit und das
% proxy-Verhalten in den Settings 1,2,4.


T = 16;
MCruns = 10;

% Beta Verteilung für lambda
alpha = 1;  %1
beta = 1;   %5

lambdaMIN = 1;
lambdaMAX = 3;

% Jacobi-Nachfrage
theta_fun = @(t) 2 + 0.5*sin(2*pi*t/10);
theta_fun = @(t) 2 + 2*(t>5);
theta_fun = @(t) 4*(0.5+0.25*sin(pi*t));    %alte Werte
d0 = 2;
d0 = 4*0.4; % alte Werte
kappa = 4;
sigma = 0.15;
dt = 0.0005;

velocityShift = 0.5;

% Inflow control
dtU = 0.5;
tiHat = 0:10;       % Unstetigkeiten von Lambda
plotBool = 0;

randeinzug = 0.1 ; % Variable für den Randeinzug damit J contain in I


% Analyse 1
dtU = 2.^(2:-1:-6);

dtU = 2.^(1:-1:-4); %vorher -7 statt -4
%dtU = [2,1,0.25,0.025]
velocityShift = [lambdaMIN,lambdaMAX];
l1 = zeros(length(dtU),1);
l2 = zeros(length(dtU),1);
l3 = zeros(length(dtU),1);
l1P = zeros(length(dtU),1);
l2P = zeros(length(dtU),1);
l3P = zeros(length(dtU),1);
l1C = zeros(length(dtU),1);
l2C = zeros(length(dtU),1);
l3C = zeros(length(dtU),1);

l1PC = zeros(length(dtU),1); % continuous proxy unten links
l2PC = zeros(length(dtU),1);
l3PC = zeros(length(dtU),1);

l1PD = zeros(length(dtU),1);  % Pwc Deterministisch oben rechts
l2PD = zeros(length(dtU),1);
l3PD = zeros(length(dtU),1);

l1D = zeros(length(dtU),1);  % continuous deterministisch oben links
l2D = zeros(length(dtU),1);
l3D = zeros(length(dtU),1);


u_all = cell(length(dtU),1);
meanOutPWC = cell(length(dtU),1);
uProxyAll = cell(length(dtU),1);
tic()
for i=1:length(dtU)
    i
    %u_all{i} = zeros(length(0:dtU(i):T),1) ;
    %u_ProxyAll{i} = zeros(length(0:dtU(i):T),1) ;
    
    % piecewise constant normal (mitte rechts)
    plotBool=1;
    [L1_error_mean,L2_error_mean,L3_error_mean,meanOutPWC{i},~,lambda,X,u_all{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,randeinzug);
    l1(i) = L1_error_mean;
    l2(i) = L2_error_mean;
    l3(i) = L3_error_mean;

    % Sample PROXY mit durchschnittlicher Geschwindigkeit und
    % durchschnittlicher Nachfrage              
    % alpha=-1 sorgt für eine deterministische Geschwindigkeit --> mit
    % Geschwindigkeit beta

    % proxy zu piecewise constant (unten rechts)
    plotBool=0;
        [L1_error_cont,L2_error_cont,L3_error_cont,mean_cont,mean_Xcont,~] = vergleicheProxy(lambda,X,dt,T,dtU(i),d0,kappa,theta_fun,alpha,beta,MCruns,i,plotBool,velocityShift,randeinzug);
        l1P(i) = L1_error_cont;
        l2P(i) = L2_error_cont;
        l3P(i) = L3_error_cont;

    %proxy ui continuous unten links
    % [L1_error_contP,L2_error_contP,L3_error_contP,mean_outP,mean_XP,lambda,XP,uP]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,1);
    %     l1PC(i) = L1_error_contP;
    %     l2PC(i) = L2_error_contP;
    %     l3PC(i) = L3_error_contP;

        % % deterministische velo pwc
        % [L1_error_contPD,L2_error_contPD,L3_error_contPD,mean_outPD,mean_XPD,lambdaPD,XPD,uPD]  = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,[2-1000*eps,2+1000*eps]);
        % l1PD(i) = L1_error_contPD;
        % l2PD(i) = L2_error_contPD;
        % l3PD(i) = L3_error_contPD;

     % Optimal ohne pwc als Referenz
     % stetige controlle mitte links
     plotBool=0;
     if i==1
        %[L1_error_cont,L2_error_cont,L3_error_cont,~,~,lambda,XCont,u_allCont{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift);
        [L1_error_cont,L2_error_cont,L3_error_cont,mean_out,mean_X,lambda,X,u]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,0,randeinzug);
        l1C(i) = L1_error_cont;
        l2C(i) = L2_error_cont;
        l3C(i) = L3_error_cont;
        l2C(1:end) = l2C(1);
        
        % stetige controlle proxy
        [L1_error_contP,L2_error_contP,L3_error_contP,mean_outP,mean_XP,lambda,XP,uP]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,1,randeinzug);
        l1PC(i) = L1_error_contP;
        l2PC(i) = L2_error_contP;
        l3PC(i) = L3_error_contP;
        l2PC(1:end) = l2PC(1);

        % % stetige controlle deterministic velocity
        % [L1_error_contD,L2_error_contD,L3_error_contD,mean_outD,mean_XD,lambdaD,X,uD]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,[2-1000*eps,2+1000*eps],1);
        % l1D(i) = L1_error_contD;
        % l2D(i) = L2_error_contD;
        % l3D(i) = L3_error_contD;
        % l2D(1:end) = l2D(1);
     end
    
end
tid = toc()

figure(2*i+1)
plot(log2(dtU), log2(l2D), LineWidth=1.5, LineStyle=":")
hold on
plot(log2(dtU), log2(l2C), LineWidth=1.5)
plot(log2(dtU), log2(l2PC), LineWidth=1.5, LineStyle="--")
plot(log2(dtU), log2(l2PD), LineWidth=1.5, LineStyle=":")
plot(log2(dtU), log2(l2), LineWidth=1.5)
plot(log2(dtU), log2(l2P), LineWidth=1.5, LineStyle="--")
legend('Continuous','Piecewise-Constant','Piecewise-Constant + Proxy')
xlabel('log2(|\Pi_n|)')
ylabel('log2(error)')
title('Logarithmic error evolution for |\Pi_n| -> 0')
legend('Continuous deterministic \lambda','Continuous','Continuous + Proxy','Piecewise-Constant deterministic \lambda','Piecewise-Constant','Piecewise-Constant + Proxy',Location='best')
hold off

figure(2*i+2)
% Errors
%y1 = log2(abs(l2D  - l2PD));   % deterministic
y2 = log2(abs(l2C  - l2));     % continuous
%y3 = log2(abs(l2PC - l2P));    % proxy
% Original curves
%h1 = plot(log2(dtU), y1, LineWidth=1.5, LineStyle=":");
hold on
h2 = plot(log2(dtU), y2, LineWidth=1.5);
h3 = plot(log2(dtU), y3, LineWidth=1.5, LineStyle="--");
% Linear fits
x = log2(dtU);
% Fit range: -3 <= log2(dtU) <= 2
idx = x >= -4 & x <= 0;
%p1 = polyfit(x(idx), y1(idx), 1);
p2 = polyfit(x(idx), y2(idx), 1);
p3 = polyfit(x(idx), y3(idx), 1);
xfit = linspace(-4, 0, 100);
% Thin dotted fits, same colors as original curves
%plot(xfit, polyval(p1,xfit), ...
%    ':', 'Color', h1.Color, 'LineWidth', 0.8)

plot(xfit, polyval(p2,xfit), ...
    ':', 'Color', h2.Color, 'LineWidth', 0.8)

plot(xfit, polyval(p3,xfit), ...
    ':', 'Color', h3.Color, 'LineWidth', 0.8)
% Labels
xlabel('$\log_2(|\Pi_n|)$', 'Interpreter', 'latex')

ylabel('$\log_2(\mathrm{error})$', 'Interpreter', 'latex')

title('Logarithmic error for decreasing the control grid cell length', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial')
% Legend
legend( ...
    '$H(u^{(n),*})-H(u^*)$', ...
    '$H(\bar{u}^{(n)})-H(\bar{u})$', ...
    ['$\mathrm{fit}: s=' sprintf('%.3f',p2(1)) '$'], ...
    ['$\mathrm{fit}: s=' sprintf('%.3f',p3(1)) '$'], ...
    'Interpreter', 'latex', ...
    'Location', 'best')
hold off

% figure(2*i+2)
% plot(log2(dtU),log2(abs(l2D - l2PD)), LineWidth=1.5, LineStyle=":")
% hold on
% plot(log2(dtU), log2(abs(l2C-l2)), LineWidth=1.5)
% plot(log2(dtU), log2(abs(l2PC - l2P)), LineWidth=1.5, LineStyle="--")
% legend('Deterministic','Stochastic','Proxy',Location='best')
% xlabel('log2(|\Pi_n|)')
% ylabel('log2(error)')
% title('Logarithmic error evolution H(u^{(n),*} - H(u^*), und determ. und proxy for |\Pi_n| -> 0')
% %legend('Continuous deterministic \lambda','Continuous','Continuous + Proxy','Piecewise-Constant deterministic \lambda','Piecewise-Constant','Piecewise-Constant + Proxy',Location='best')
% hold off


% figure(2*i+2)
% hold on
% for j=1:i
%     plot(0:dtU(j):T-dtU(j), uProxyAll{j}-u_all{j}, 'LineWidth',1.5)
% end
% hold off
% xlabel('time')
% ylabel('inflow difference')
% labels = compose("\\Delta t_u = %g", dtU);
% legend(labels)
% title('Difference of Proxy inflow and MC mean inflow of optimal solution')

% plot ui
%dtU = 2.^(2:-1:-5);
leg_entries = cell(1,length(1:2:length(u_all))+1);
figure(37)
for k = 1:2:length(u_all)
hold on
    ui = u_all{k};
    n = length(ui);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, ui, 'LineWidth', 1.2);
    leg_entries{(k+1)/2} = ['k = ' num2str(k)];
hold off
end

hold on
plot(dt:dt:T,u,'LineWidth',2)
leg_entries{end} = 'continuous';
%plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'Marker','*')


hold off
xlabel('t');
ylabel('u');
title('Comparison of optimal inflows');
legend(leg_entries, 'Location', 'best');
grid on;

figure(38)
for k = 1:2:length(meanOutPWC)
hold on
    outi = meanOutPWC{k};
    n = length(outi);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, outi, 'LineWidth', 1.2);
hold off
end

hold on
plot(0:dt:T,mean_out,'LineWidth',2,'LineStyle',':')
plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'LineStyle','--')


hold off
xlabel('t');
ylabel('out');
title('Comparison outflows and demand');
leg_entries{end+1} = 'mean demand';
legend(leg_entries, 'Location', 'best');
grid on;


%% Decreasing the variance of lambda (Figure 4,5)
% aber variablem lambdaMIN und lambdaMAX aber zentriert in 2

T = 16;
MCruns = 10;
%print('Lauf 2')
% Beta Verteilung für lambda
alpha = 1;  %1
beta = 1;   %5

%lambdaMIN = 1;
%lambdaMAX = 3;

% Jacobi-Nachfrage
theta_fun = @(t) 2 + 0.5*sin(2*pi*t/10);
theta_fun = @(t) 2 + 2*(t>5);
theta_fun = @(t) 4*(0.5+0.25*sin(pi*t));    %alte Werte
d0 = 2;
d0 = 4*0.4; % alte Werte
kappa = 4;
sigma = 0.15;
dt = 0.0005;

velocityShift = 0.5;
randeinzug = 2 ; % Variable für den Randeinzug damit J contain in I


% Inflow control
dtU = 0.5;
tiHat = 0:10;       % Unstetigkeiten von Lambda
plotBool = 0;


% Analyse 2

shifts = (3:-1:-6);

lambdaMIN = 2-0.5*sqrt(2.^shifts);
lambdaMAX = 2+0.5*sqrt(2.^shifts);
%dtU = [2,1,0.25,0.025]

l1 = zeros(length(shifts),1); % pwc mitte rechts
l2 = zeros(length(shifts),1);
l3 = zeros(length(shifts),1);

l1P = zeros(length(shifts),1); % Proxy pwc unten rechts
l2P = zeros(length(shifts),1);
l3P = zeros(length(shifts),1);

l1C = zeros(length(shifts),1); % continuous mitte links
l2C = zeros(length(shifts),1);
l3C = zeros(length(shifts),1);

l1PC = zeros(length(shifts),1); % continuous proxy unten links
l2PC = zeros(length(shifts),1);
l3PC = zeros(length(shifts),1);

l1PD = zeros(length(shifts),1);  % Pwc Deterministisch oben rechts
l2PD = zeros(length(shifts),1);
l3PD = zeros(length(shifts),1);

l1D = zeros(length(shifts),1);  % continuous deterministisch oben links
l2D = zeros(length(shifts),1);
l3D = zeros(length(shifts),1);


u_all = cell(length(shifts),1);

u_allPWC = cell(length(shifts),1);
meanOutPWC = cell(length(shifts),1);
meanOutC = cell(length(shifts),1);
uProxyAll = cell(length(shifts),1);
tic()
for i=1:length(shifts)
    i
    %u_all{i} = zeros(length(0:dtU(i):T),1) ;
    %u_ProxyAll{i} = zeros(length(0:dtU(i):T),1) ;

    %Stochastic Piecewise constant
    velocityShift = [lambdaMIN(i),lambdaMAX(i)];
    plotBool=0;
    [L1_error_mean,L2_error_mean,L3_error_mean,meanOutPWC{i},~,lambda,X,u_allPWC{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,i,velocityShift,randeinzug);
    l1(i) = L1_error_mean;
    l2(i) = L2_error_mean;
    l3(i) = L3_error_mean;

    % Sample PROXY mit durchschnittlicher Geschwindigkeit und
    % durchschnittlicher Nachfrage              
    % 
plotBool=0;
    [L1_error_cont,L2_error_cont,L3_error_cont,mean_cont,mean_Xcont,uPP] = vergleicheProxy(lambda,X,dt,T,dtU,d0,kappa,theta_fun,alpha,beta,MCruns,i,plotBool,velocityShift,randeinzug);
    l1P(i) = L1_error_cont;
    l2P(i) = L2_error_cont;
    l3P(i) = L3_error_cont;

     % Optimal ohne pwc als Referenz
     plotBool=0;
    % if i==1
        %[L1_error_cont,L2_error_cont,L3_error_cont,~,~,lambda,XCont,u_allCont{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift);
        [L1_error_cont,L2_error_cont,L3_error_cont,meanOutC{i},mean_X,lambda,X,u_all{i}]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,i,velocityShift,0,randeinzug);
        l1C(i) = L1_error_cont;
        l2C(i) = L2_error_cont;
        l3C(i) = L3_error_cont;
     %   l2C(1:end) = l2C(1);
     %end
    [L1_error_contP,L2_error_contP,L3_error_contP,mean_outP,mean_XP,lambda,XP,uP]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,i,velocityShift,1,randeinzug);
        l1PC(i) = L1_error_contP;
        l2PC(i) = L2_error_contP;
        l3PC(i) = L3_error_contP;
     % Optimal ohne pwc und DETERMINISTISCHER Velo als Referenz
     plotBool=0;
     if i==1
        % %[L1_error_cont,L2_error_cont,L3_error_cont,~,~,lambda,XCont,u_allCont{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift);
        % [L1_error_contD,L2_error_contD,L3_error_contD,mean_outD,mean_XD,lambdaD,X,uD]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,i,[2-1000*eps,2+1000*eps],0,randeinzug);
        % l1D(i) = L1_error_contD;
        % l2D(i) = L2_error_contD;
        % l3D(i) = L3_error_contD;
        % l2D(1:end) = l2D(1);

        % [L1_error_contPD,L2_error_contPD,L3_error_contPD,mean_outPD,mean_XPD,lambdaPD,XP,uPD]  = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,i,[2-1000*eps,2+1000*eps],randeinzug);
        % l1PD(i) = L1_error_contPD;
        % l2PD(i) = L2_error_contPD;
        % l3PD(i) = L3_error_contPD;
        % l2PD(1:end) = l2PD(1);
     end
    
end
tid = toc()

figure(49)
%plot(log2(dtU), log2(l1), LineWidth=1.5)
%plot((shifts), log2(l2D), LineWidth=1.5)
hold on
plot((shifts), log2(l2C), LineWidth=1.5)
plot((shifts), log2(l2PC), LineWidth=1.5)
%plot((shifts), log2(l2PD), LineWidth=1.5, LineStyle="--")
plot((shifts), log2(l2), LineWidth=1.5,LineStyle="--")
%plot(log2(dtU), log2(l3), LineWidth=1.5)
%plot(log2(dtU), log2(l1P), LineWidth=1.5, LineStyle="--")

plot((shifts), log2(l2P), LineWidth=1.5, LineStyle="--")

%plot(log2(dtU), log2(l3P), LineWidth=1.5, LineStyle="--")
%legend('L1-error','L2-error','L3-error','L1-errorProxy','L2-errorProxy','L3-errorProxy')
legend('Continuous deterministic \lambda','Continuous','Continuous + Proxy','Piecewise-Constant deterministic \lambda','Piecewise-Constant','Piecewise-Constant + Proxy',Location='best')
xlabel('log2(Var(\lambda))')
ylabel('log2(error)')
title('Logarithmic error evolution for Var(\lambda) -> 0')
hold off


figure(50)
% Data
y1 = log2(abs(l2PC - l2C));
y2 = log2(abs(l2P  - l2));
% Original curves
h1 = plot(shifts, y1, LineWidth=2);
hold on
h2 = plot(shifts, y2, LineWidth=2, LineStyle="--");
% Linear fits for shifts = -2,...,4
idx = shifts >= -3 & shifts <= 3;
p1 = polyfit(shifts(idx), y1(idx), 1);
p2 = polyfit(shifts(idx), y2(idx), 1);
xfit = linspace(-3, 3, 100);
% Fits with same colors as original curves
plot(xfit, polyval(p1, xfit), ...
    ':', 'Color', h1.Color, 'LineWidth', 1.3)
plot(xfit, polyval(p2, xfit), ...
    ':', 'Color', h2.Color, 'LineWidth', 1.3)
xlabel('$\log_2(\mathrm{Var}(\lambda))$', 'Interpreter', 'latex')
ylabel('$\log_2(\mathrm{error})$', 'Interpreter', 'latex')
title('Logarithmic error for decreasing variance in $\lambda$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial')
legend( ...
    '$\log(H(\bar{u})-H(u^*))$', ...
    '$\log(H(\bar{u}^{(n)})-H(u^{(n),*}))~~$', ...
    ['$\mathrm{fit}: slope=' sprintf('%.3f',p1(1)) '$'], ...
    ['$\mathrm{fit}: slope=' sprintf('%.3f',p2(1)) '$'], ...
    'Interpreter', 'latex', ...
    'Location', 'best')
hold off

figure(51)
% Data
y1 = log2(abs(l2PC - l2C));
y2 = log2(abs(l2P  - l2));
% Original curves
%h1 = plot(shifts, y1, LineWidth=2);
h1 = plot(shifts(idx), y1(idx), LineWidth=2);
hold on
%h2 = plot(shifts, y2, LineWidth=2, LineStyle="--");
h2 = plot(shifts(idx), y2(idx), LineWidth=2, LineStyle="--");
% Linear fits for shifts = -2,...,4
idx = shifts >= -2 & shifts <= 4;
p1 = polyfit(shifts(idx), y1(idx), 1);
p2 = polyfit(shifts(idx), y2(idx), 1);
xfit = linspace(-2, 4, 100);
% Fits with same colors as original curves
plot(xfit, polyval(p1, xfit), ...
    ':', 'Color', h1.Color, 'LineWidth', 1.3)
plot(xfit, polyval(p2, xfit), ...
    ':', 'Color', h2.Color, 'LineWidth', 1.3)
xlabel('$\log_2(\mathrm{Var}(\lambda))$', 'Interpreter', 'latex')
ylabel('$\log_2(\mathrm{error})$', 'Interpreter', 'latex')
title('Logarithmic error for decreasing variance in $\lambda$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial')
legend( ...
    '$\log(H(\bar{u})-H(u^*))$', ...
    '$\log(H(\bar{u}^{(n)})-H(u^{(n),*}))~~$', ...
    ['$\mathrm{fit}: slope=' sprintf('%.3f',p1(1)) '$'], ...
    ['$\mathrm{fit}: slope=' sprintf('%.3f',p2(1)) '$'], ...
    'Interpreter', 'latex', ...
    'Location', 'best')
hold off





% plot ui
%shifts = (4:-1:-8);

lambdaMIN = 2-0.5*sqrt(2.^shifts);
lambdaMAX = 2+0.5*sqrt(2.^shifts);
leg_entries = cell(1,length(1:2:length(u_all))+2);
figure(73)
for k = 1:2:length(u_all)
hold on
    ui = u_all{k};
    n = length(ui);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, ui, 'LineWidth', 1.2);
    leg_entries{(k+1)/2} = ['k = ' num2str(shifts(k))];
hold off
end

hold on
plot(dt:dt:T,uD,'LineWidth',2, 'LineStyle','--')
leg_entries{end-1} = 'deterministic';
%plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'Marker','*')
plot(dt:dt:T,uP,'LineWidth',2,'LineStyle',':')
leg_entries{end} = 'proxy';

hold off
xlabel('t');
ylabel('u');
title('Comparison of the optimal inflows for continuous controls');
legend(leg_entries, 'Location', 'best');
grid on;



figure(74)
for k = 1:2:length(meanOutC)
hold on
    outi = meanOutC{k};
    n = length(outi);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, outi, 'LineWidth', 1.2);
hold off
end

hold on
plot(0:dt:T,mean_outD,'LineWidth',2,'LineStyle',':')
plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'LineStyle','--')


hold off
xlabel('t');
ylabel('outdlow');
title('Comparison of the utflows and demand for continuous controls');
leg_entries{end} = 'mean demand';
legend(leg_entries, 'Location', 'best');
grid on;

print('Hier einmal prüfen ob die richtigen Kurven geplottet werden!')

% nochmal piecewise constant


lambdaMIN = 2-0.5*sqrt(2.^shifts);
lambdaMAX = 2+0.5*sqrt(2.^shifts);
leg_entries = cell(1,length(1:2:length(u_allPWC))+2);


figure(94)
for k = 1:2:length(u_allPWC)
hold on
    ui = u_allPWC{k};
    n = length(ui);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(0:dtU:(T-dtU), ui, 'LineWidth', 1.2);
    leg_entries{(k+1)/2} = ['k = ' num2str(shifts(k))];
hold off
end

hold on
stairs(0:dtU:(T-dtU),uPD,'LineWidth',2, 'LineStyle','--')
leg_entries{end-1} = 'deterministic';
%plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'Marker','*')
stairs(0:dtU:(T-dtU),uPP,'LineWidth',2,'LineStyle',':')
leg_entries{end} = 'proxy';

hold off
xlabel('t');
ylabel('u');
title('Comparison of the optimal inflows for piecewise constant control');
legend(leg_entries, 'Location', 'best');
grid on;



figure(95)
for k = 1:2:length(meanOutPWC)
hold on
    outi = meanOutPWC{k};
    n = length(outi);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, outi, 'LineWidth', 1.2);
hold off
end

hold on
plot(0:dt:T,mean_outPD,'LineWidth',2,'LineStyle',':')
plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'LineStyle','--')


hold off
xlabel('t');
ylabel('outflow');
title('Comparison of the outflows and demand for piecewise constant control');
leg_entries{end} = 'mean demand';
legend(leg_entries, 'Location', 'best');
grid on;



%% Boundary effects (Figure 6, left)


T = 6;
MCruns = 25000;

lambdaMIN = 1;
lambdaMAX = 3;
velocityShift = [lambdaMIN,lambdaMAX]
tic()
dtU = 0.25;
proxy = 0;
randeinzug = 0;
ii=1;
[L1_error_mean,L2_error_mean,L3_error_mean,mean_out,mean_X,lambda,X,u]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,ii,velocityShift,proxy,randeinzug);

[L1_error_meanWB,L2_error_meanWB,L3_error_meanWB,mean_outWB,mean_XWB,lambdaWB,XWB,uWB]  = sampleSingleSpeedMCcontinuous_wrongBoundary(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU,plotBool,ii,velocityShift,proxy,randeinzug);


tid = toc()

figure(154)

plot(0:dt:(T-(1/3)),u(1:round(5.66666/dt)+1),'LineWidth',2)
hold on
plot(0:dt:(T-1/3),uWB(1:round(5.66666/dt)+1),'LineWidth',2,'LineStyle','--')
legend('optimal','interior')
xlim([0 5.6666666])
xlabel('time')
ylabel('controls')
title('Comparison of optimal control and unconditioned interior control')
hold off

%% Boundary effects (Figure 6, right)
T = 6;
MCruns = 10;

% Beta Verteilung für lambda
alpha = 1;  %1
beta = 1;   %5

lambdaMIN = 1;
lambdaMAX = 3;

% Jacobi-Nachfrage
theta_fun = @(t) 2 + 0.5*sin(2*pi*t/10);
theta_fun = @(t) 2 + 2*(t>5);
theta_fun = @(t) 4*(0.5+0.25*sin(pi*t));    %alte Werte
d0 = 2;
d0 = 4*0.4; % alte Werte
kappa = 4;
sigma = 0.15;
dt = 0.0005;

velocityShift = 0.5;

% Inflow control
dtU = 0.5;
tiHat = 0:10;       % Unstetigkeiten von Lambda
plotBool = 0;

randeinzug = 0.1 ; % Variable für den Randeinzug damit J contain in I


% Analyse 1
dtU = 2.^(2:-1:-6);

dtU = 2.^(1:-1:-4); %vorher -7 statt -4
%dtU = [2,1,0.25,0.025]
velocityShift = [lambdaMIN,lambdaMAX];
l1 = zeros(length(dtU),1);
l2 = zeros(length(dtU),1);
l3 = zeros(length(dtU),1);
l1P = zeros(length(dtU),1);
l2P = zeros(length(dtU),1);
l3P = zeros(length(dtU),1);
l1C = zeros(length(dtU),1);
l2C = zeros(length(dtU),1);
l3C = zeros(length(dtU),1);

l1PC = zeros(length(dtU),1); % continuous proxy unten links
l2PC = zeros(length(dtU),1);
l3PC = zeros(length(dtU),1);

l1PD = zeros(length(dtU),1);  % Pwc Deterministisch oben rechts
l2PD = zeros(length(dtU),1);
l3PD = zeros(length(dtU),1);

l1D = zeros(length(dtU),1);  % continuous deterministisch oben links
l2D = zeros(length(dtU),1);
l3D = zeros(length(dtU),1);


u_all = cell(length(dtU),1);
meanOutPWC = cell(length(dtU),1);
uProxyAll = cell(length(dtU),1);
tic()
for i=1:length(dtU)
    i
    %u_all{i} = zeros(length(0:dtU(i):T),1) ;
    %u_ProxyAll{i} = zeros(length(0:dtU(i):T),1) ;
    
    % piecewise constant normal (mitte rechts)
    plotBool=1;
    [L1_error_mean,L2_error_mean,L3_error_mean,meanOutPWC{i},~,lambda,X,u_all{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,randeinzug);
    l1(i) = L1_error_mean;
    l2(i) = L2_error_mean;
    l3(i) = L3_error_mean;

    % Sample PROXY mit durchschnittlicher Geschwindigkeit und
    % durchschnittlicher Nachfrage              
    % alpha=-1 sorgt für eine deterministische Geschwindigkeit --> mit
    % Geschwindigkeit beta

    % proxy zu piecewise constant (unten rechts)
    plotBool=0;
        [L1_error_cont,L2_error_cont,L3_error_cont,mean_cont,mean_Xcont,~] = vergleicheProxy(lambda,X,dt,T,dtU(i),d0,kappa,theta_fun,alpha,beta,MCruns,i,plotBool,velocityShift,randeinzug);
        l1P(i) = L1_error_cont;
        l2P(i) = L2_error_cont;
        l3P(i) = L3_error_cont;

    %proxy ui continuous unten links
    % [L1_error_contP,L2_error_contP,L3_error_contP,mean_outP,mean_XP,lambda,XP,uP]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,1);
    %     l1PC(i) = L1_error_contP;
    %     l2PC(i) = L2_error_contP;
    %     l3PC(i) = L3_error_contP;

        % % deterministische velo pwc
        % [L1_error_contPD,L2_error_contPD,L3_error_contPD,mean_outPD,mean_XPD,lambdaPD,XPD,uPD]  = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,[2-1000*eps,2+1000*eps]);
        % l1PD(i) = L1_error_contPD;
        % l2PD(i) = L2_error_contPD;
        % l3PD(i) = L3_error_contPD;

     % Optimal ohne pwc als Referenz
     % stetige controlle mitte links
     plotBool=0;
     if i==1
        %[L1_error_cont,L2_error_cont,L3_error_cont,~,~,lambda,XCont,u_allCont{i},~,~] = sampleSingleSpeedMC(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift);
        [L1_error_cont,L2_error_cont,L3_error_cont,mean_out,mean_X,lambda,X,u]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,0,randeinzug);
        l1C(i) = L1_error_cont;
        l2C(i) = L2_error_cont;
        l3C(i) = L3_error_cont;
        l2C(1:end) = l2C(1);
        
        % stetige controlle proxy
        [L1_error_contP,L2_error_contP,L3_error_contP,mean_outP,mean_XP,lambda,XP,uP]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,velocityShift,1,randeinzug);
        l1PC(i) = L1_error_contP;
        l2PC(i) = L2_error_contP;
        l3PC(i) = L3_error_contP;
        l2PC(1:end) = l2PC(1);

        % % stetige controlle deterministic velocity
        % [L1_error_contD,L2_error_contD,L3_error_contD,mean_outD,mean_XD,lambdaD,X,uD]  = sampleSingleSpeedMCcontinuous(T,alpha,beta,MCruns,theta_fun,d0,kappa,sigma,dt,dtU(i),plotBool,i,[2-1000*eps,2+1000*eps],1);
        % l1D(i) = L1_error_contD;
        % l2D(i) = L2_error_contD;
        % l3D(i) = L3_error_contD;
        % l2D(1:end) = l2D(1);
     end
    
end
tid = toc()

figure(2*i+1)
plot(log2(dtU), log2(l2D), LineWidth=1.5, LineStyle=":")
hold on
plot(log2(dtU), log2(l2C), LineWidth=1.5)
plot(log2(dtU), log2(l2PC), LineWidth=1.5, LineStyle="--")
plot(log2(dtU), log2(l2PD), LineWidth=1.5, LineStyle=":")
plot(log2(dtU), log2(l2), LineWidth=1.5)
plot(log2(dtU), log2(l2P), LineWidth=1.5, LineStyle="--")
legend('Continuous','Piecewise-Constant','Piecewise-Constant + Proxy')
xlabel('log2(|\Pi_n|)')
ylabel('log2(error)')
title('Logarithmic error evolution for |\Pi_n| -> 0')
legend('Continuous deterministic \lambda','Continuous','Continuous + Proxy','Piecewise-Constant deterministic \lambda','Piecewise-Constant','Piecewise-Constant + Proxy',Location='best')
hold off

figure(2*i+2)
% Errors
%y1 = log2(abs(l2D  - l2PD));   % deterministic
y2 = log2(abs(l2C  - l2));     % continuous
%y3 = log2(abs(l2PC - l2P));    % proxy
% Original curves
%h1 = plot(log2(dtU), y1, LineWidth=1.5, LineStyle=":");
hold on
h2 = plot(log2(dtU), y2, LineWidth=1.5);
h3 = plot(log2(dtU), y3, LineWidth=1.5, LineStyle="--");
% Linear fits
x = log2(dtU);
% Fit range: -3 <= log2(dtU) <= 2
idx = x >= -4 & x <= 0;
%p1 = polyfit(x(idx), y1(idx), 1);
p2 = polyfit(x(idx), y2(idx), 1);
p3 = polyfit(x(idx), y3(idx), 1);
xfit = linspace(-4, 0, 100);
% Thin dotted fits, same colors as original curves
%plot(xfit, polyval(p1,xfit), ...
%    ':', 'Color', h1.Color, 'LineWidth', 0.8)

plot(xfit, polyval(p2,xfit), ...
    ':', 'Color', h2.Color, 'LineWidth', 0.8)

plot(xfit, polyval(p3,xfit), ...
    ':', 'Color', h3.Color, 'LineWidth', 0.8)
% Labels
xlabel('$\log_2(|\Pi_n|)$', 'Interpreter', 'latex')

ylabel('$\log_2(\mathrm{error})$', 'Interpreter', 'latex')

title('Logarithmic error for decreasing the control grid cell length', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial')
% Legend
legend( ...
    '$H(u^{(n),*})-H(u^*)$', ...
    '$H(\bar{u}^{(n)})-H(\bar{u})$', ...
    ['$\mathrm{fit}: s=' sprintf('%.3f',p2(1)) '$'], ...
    ['$\mathrm{fit}: s=' sprintf('%.3f',p3(1)) '$'], ...
    'Interpreter', 'latex', ...
    'Location', 'best')
hold off

% figure(2*i+2)
% plot(log2(dtU),log2(abs(l2D - l2PD)), LineWidth=1.5, LineStyle=":")
% hold on
% plot(log2(dtU), log2(abs(l2C-l2)), LineWidth=1.5)
% plot(log2(dtU), log2(abs(l2PC - l2P)), LineWidth=1.5, LineStyle="--")
% legend('Deterministic','Stochastic','Proxy',Location='best')
% xlabel('log2(|\Pi_n|)')
% ylabel('log2(error)')
% title('Logarithmic error evolution H(u^{(n),*} - H(u^*), und determ. und proxy for |\Pi_n| -> 0')
% %legend('Continuous deterministic \lambda','Continuous','Continuous + Proxy','Piecewise-Constant deterministic \lambda','Piecewise-Constant','Piecewise-Constant + Proxy',Location='best')
% hold off


% figure(2*i+2)
% hold on
% for j=1:i
%     plot(0:dtU(j):T-dtU(j), uProxyAll{j}-u_all{j}, 'LineWidth',1.5)
% end
% hold off
% xlabel('time')
% ylabel('inflow difference')
% labels = compose("\\Delta t_u = %g", dtU);
% legend(labels)
% title('Difference of Proxy inflow and MC mean inflow of optimal solution')

% plot ui
%dtU = 2.^(2:-1:-5);
leg_entries = cell(1,length(1:2:length(u_all))+1);
figure(37)
for k = 1:2:length(u_all)
hold on
    ui = u_all{k};
    n = length(ui);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, ui, 'LineWidth', 1.2);
    leg_entries{(k+1)/2} = ['k = ' num2str(k)];
hold off
end

hold on
plot(dt:dt:T,u,'LineWidth',2)
leg_entries{end} = 'continuous';
%plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'Marker','*')


hold off
xlabel('t');
ylabel('u');
title('Comparison of optimal inflows');
legend(leg_entries, 'Location', 'best');
grid on;

figure(38)
for k = 1:2:length(meanOutPWC)
hold on
    outi = meanOutPWC{k};
    n = length(outi);

    % äquidistante Zeitpunkte
    t = linspace(0, T, n);

    % stairs-Plot
    stairs(t, outi, 'LineWidth', 1.2);
hold off
end

hold on
plot(0:dt:T,mean_out,'LineWidth',2,'LineStyle',':')
plot(dt:dt:T, jacobiMeanClosed(dt:dt:T, d0, kappa, theta_fun),'LineWidth',2,'LineStyle','--')


hold off
xlabel('t');
ylabel('out');
title('Comparison outflows and demand');
leg_entries{end+1} = 'mean demand';
legend(leg_entries, 'Location', 'best');
grid on;

