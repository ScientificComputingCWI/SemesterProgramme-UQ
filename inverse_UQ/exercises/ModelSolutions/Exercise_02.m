clear all; clc; close all;

%% (I) Forward map and its discretization
%
r = 0.25;
x0 = 0.1;

% Time of observation and their corresponding indices 
t_obs = 1:3;

% RHS of ODE
f_k = @(t,x,k) r * x * (k - x); 

% Forward map: Euler discretization of ODE and observation map
Gh = @(k,h) Euler(h, x0, 3, @(t,x) f_k(t,x,k), t_obs);

% True forward map without discretization
G = @(k) k ./ (1 + exp(-r*k* t_obs) * (k/x0-1));

%% Data and true posterior
% Discretizing parameter domain
hk = 0.01;
ks = 10:hk:20;

% Observational data
y = [3, 14, 17];

% Noise variance
sigma2 = 1;

% True Potential
Pot = @(k) norm(y - G(k)).^2/2/sigma2;

% Prior density
pi0_fun = @(k) 1/10;

% Computing posterior density at grid points ks
piy = zeros(length(ks),1);
tic;
for i = 1:length(ks),
    piy(i) = exp(- Pot(ks(i))) * pi0_fun(ks(i));
end
toc;

% Computing normalizing constant via trapezoidal rule
% quadrature weights
w = [0.5, ones(1, length(ks)-2), 0.5] / hk;
Z = w * piy;
piy = piy/Z;

figure(1);
plot(ks, piy,'-k'); 
title('Posteriors')
hold on;
grid on;
names{1} = 'true'
set(gca,'Fontsize',14)

%% Approximation of posterior
hs = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];

% Approximative Potential
Poth = @(k,h) norm(y - Gh(k,h)).^2/2/sigma2;

dHell = zeros(1, length(hs));

for j = 1:length(hs),
    piyh = zeros(length(ks),1);
    tic;
    for i = 1:length(ks),
        piyh(i) = exp(- Poth(ks(i), hs(j))) * pi0_fun(ks(i));
    end
    Z = w * piyh;
    piyh = piyh/Z;
    
    % Plotting
    if mod(j,2) == 1,
        figure(1);
        plot(ks, piyh,'--'); 
    end

    dHell(j) = sqrt( w * ( sqrt(piy) - sqrt(piyh)).^2 );
end

% Adding legend to figure(1)
names{2} = sprintf('h = %.3f', hs(1));
names{3} = sprintf('h = %.3f', hs(3));
names{4} = sprintf('h = %.3f', hs(5));
names{5} = sprintf('h = %.3f', hs(7));
names{6} = sprintf('h = %.3f', hs(9));
figure(1)
legend(names)
set(gca,'fontsize',14)

% Plotting decay in Hellinger distance
figure();
loglog(hs, dHell,'-ob')
title('Error in Hellinger distance')
xlabel('h')
ylabel('Helinger distance')
grid on
set(gca,'Fontsize',14)

return
%% (I.1)  Illustrations for numerical approximation
k_fun = @(t) 15 * exp(1/8 * sqrt(2) * (1.5*sin(pi*t) + 0.5*sin(2*pi*t)) );
k = k_fun(ts);
N = F(k);
y = ObsOp(k);

figure();
plot(ts, k,'-k','LineWidth',2);
grid on;
xlabel('t')
ylabel('k(t)')
ylim([0,22])
set(gca,'fontsize',14)

figure();
plot(ts, N,'-r','LineWidth',2);
grid on;
xlabel('t')
ylabel('N(t)')
ylim([0,22])
set(gca,'fontsize',14)

figure();
plot(t_obs, y, 'ob','LineWidth',2);
grid on;
xlabel('t')
ylim([0,22])
set(gca,'fontsize',14)

figure();
plot(ts(1:50:end), N(1:50:end), 'o--r','LineWidth',2,'MarkerFaceColor','r');
grid on;
xlabel('t')
ylabel('N(t)')
ylim([0,22])
set(gca,'fontsize',14)



%% (II) Bayesian inverse problem
%

%% (II.1) Ground truth and data
k_true = 15 * exp(1/8 * sqrt(2)*sin(2*pi*ts));

N_true = F(k_true);
y_true = ObsOp(k_true);

sigma = 0.25;
y_obs = y_true + sigma * randn(length(y_true),1);

%% Plotting

% Function k
figure();
plot(ts, k_true,'-k','LineWidth',2);
grid on;
xlabel('t')
ylim([0,20])
legend('Kapazitaet k^*','Location','SouthEast');
set(gca,'fontsize',14)

% Function k and solution N 
figure();
plot(ts, k_true,'-k','LineWidth',2);
hold on;
plot(ts, N_true,'-r','LineWidth',2);
grid on;
xlabel('t')
ylim([0,20])
legend('Kapazitaet k^*','Population N^*','Location','SouthEast');
set(gca,'fontsize',14)

% Function k, solution N  and observations
figure();
plot(ts, k_true,'-k','LineWidth',2);
hold on;
plot(ts, N_true,'-r','LineWidth',2);
plot(t_obs, y_true,'ob','LineWidth',2);
grid on;
xlabel('t')
ylim([0,20])
legend('Kapazitaet k^*','Population N^*','Beobachtungen y^*','Location','SouthEast');
set(gca,'fontsize',14)

% Function k, solution N and noisy observations
figure();
plot(ts, k_true,'-k','LineWidth',2);
hold on;
plot(ts, N_true,'-r','LineWidth',2);
plot(t_obs, y_obs,'ob','LineWidth',2);
grid on;
xlabel('t')
ylim([0,20])
legend('Kapazitaet k^*','Population N^*','Beobachtungen y','Location','SouthEast');
set(gca,'fontsize',14)

% Function k and noisy observations
figure();
plot(ts, k_true,'--','LineWidth',2,'Color',[.75 .75 .75]);
grid on;
xlabel('t')
ylim([0,20])
hold on;
plot(t_obs, y_obs,'ob','LineWidth',2);
set(gca,'fontsize',14)

%% (II.2) Prior

M_KL = 500;

% Decay of KL eigenvalues
p = 5; % p = 2, 3, 4 oder 5
lambda = 1 ./ (pi * (1:M_KL)).^p;

% Normalization for totatl variance
tot_var = 0.01;
lambda = tot_var/sum(lambda) * lambda;

% KL eigenfunctions
KL = zeros(length(ts),M_KL);
for i = 1:M_KL,
    KL(:,i) = sqrt(2) * sin(i * pi * ts');
end

% Mean such  that lognormal field has specific mean 
mean_KL = -sqrt(2)/2 * ( KL.^2 * lambda');

% Generation of k via KL system
k_KL = @(xi) 15 * exp(mean_KL + KL * xi);

% Prior covariance
C0 = diag(lambda);
L0 = diag(sqrt(lambda));

%% Plotting prior

% 100 realizations and truth
figure();
plot(ts, k_KL(L0 * randn(M_KL,100)),'LineWidth',1); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title({'100 Realisierungen des Priors','und wahre Kapazitaet (schwarz)'})
set(gca,'Fontsize',14)

y_limits = get(gca,'ylim');

% 5 realizations and truth
figure();
plot(ts, k_KL(L0 * randn(M_KL,5)),'LineWidth',2); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title({'5 Realisierungen des Priors','und wahre Kapazitaet (schwarz)'})
set(gca,'Fontsize',14)
ylim(y_limits)

% Quantile-Plot
Ks = k_KL(L0 * randn(M_KL,1e5));
quant_005 = quantile(Ks,0.05,2);
quant_095 = quantile(Ks,0.95,2);
mean_Ks = mean(Ks,2);
clear Ks;
figure();
plot(ts, mean_Ks,'-b','LineWidth',2); hold on
plot(ts, quant_005,'--r','LineWidth',2);
plot(ts, quant_095,'--r','LineWidth',2);
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
legend('Mittelwert','5%-Quantil','95%-Quantil','wahre Kapazitaet','Location','NorthEast')
title({'A-priori Mittelwert und Quantile','und wahre Kapazitaet (schwarz)'})
set(gca,'Fontsize',14)

%% (II.3) MCMC - pCN

% Potential operator
PotOp = @(xi) norm(y_obs - ObsOp(k_KL(xi)))^2 / sigma^2;

M_pCN = 1e5; %lenght of chain
Xis = zeros(M_KL, M_pCN); % states of chain

% Stepsize
s = 0.15;

% Pre-tuned values:
% s = 0.15; % p = 2;
% s = 0.15; % p = 3;
% s = 0.15; % p = 4;
% s = 0.15; % p = 5;

%% Burn-in phase
M_burnin = round(0.1 * M_pCN);
xi = zeros(M_KL,1);
xi_pot = PotOp(xi);
A = 0;
tic;
for i = 1:M_burnin,
    eta = sqrt(1-s^2) * xi + s * L0 * randn(M_KL,1);
    eta_pot = PotOp(eta);
    
    a = min(1, exp(xi_pot - eta_pot));
    A = A + a;
    
    if rand(1) <= a,
        xi = eta;
        xi_pot = eta_pot;
    end
end
toc;
A/M_burnin

%% Real run
A = 0;
tic;
for i = 1:M_pCN,
    eta = sqrt(1-s^2) * xi + s * L0 * randn(M_KL,1);
    eta_pot = PotOp(eta);
    
    a = min(1, exp(xi_pot - eta_pot));
    A = A + a;
    
    if rand(1) <= a,
        xi = eta;
        xi_pot = eta_pot;        
    end
    
    Xis(:,i) = xi;
end
toc;
A/M_pCN

%% Plotting posterior

% 100 realizations and truth
figure();
plot(ts, k_KL(Xis(:,1:(M_pCN/100):end) ),'LineWidth',1); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title({'100 Realisierungen des Posteriors','und wahre Kapazitaet (schwarz)'})
set(gca,'Fontsize',14)

y_limits = get(gca,'ylim');

% 5 realizations and truth
figure();
plot(ts, k_KL( Xis(:,1:(M_pCN/4):end) ),'LineWidth',2); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title({'5 Realisierungen des Posteriors','und wahre Kapazitaet (schwarz)'})
set(gca,'Fontsize',14)
ylim(y_limits)

% Quantile-Plot
Ks = k_KL(Xis);
quant_005 = quantile(Ks,0.05,2);
quant_095 = quantile(Ks,0.95,2);
mean_Ks = mean(Ks,2);
clear Ks;
figure();
plot(ts, mean_Ks,'-b','LineWidth',2); hold on
plot(ts, quant_005,'--r','LineWidth',2);
plot(ts, quant_095,'--r','LineWidth',2);
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
legend('Mittelwert','5%-Quantil','95%-Quantil','wahre Kapazitaet','Location','NorthEast')
title({'A-posteriori Mittelwert und Quantile','und wahre Kapazitaet (schwarz)'})
set(gca,'Fontsize',14)

% Plot Path of Markov chain
figure();
plot(1:2000,Xis(1:3,1:2000),'.-','LineWidth',2);
grid on;
xlabel('Iteration der Markovkette')
legend('\xi_1','\xi_2','\xi_3')
set(gca,'FontSize',14)
title('Pfad der Markovkette')


% Compare prior and posterior for certain xi
for j = [1,2,3,10,100],
    figure();
    histogram(Xis(j,:),'normalization','pdf'); hold on;
    xi_grid = -3*sqrt(lambda(j)):6*sqrt(lambda(j))/1000:3*sqrt(lambda(j));
    plot(xi_grid,normpdf(xi_grid,0,sqrt(lambda(j))),'-r','LineWidth',2)
    xlabel(['\xi_',sprintf('{%i}',j)])
    grid on;
    legend('a-posteriori','a-priori')
    set(gca,'FontSize',14)
end