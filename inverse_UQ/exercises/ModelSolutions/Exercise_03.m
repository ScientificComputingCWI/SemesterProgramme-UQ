%clear all; clc; close all;

%% (I) Forward map and its discretization
%
r = 0.25;
x0 = 1;

% Time of observation and their corresponding indices 
t_obs = 1:12;

% time discretization
h = 0.01;
ts = 0:h:12;

% Forward map: Euler discretization of ODE and observation map
Gh = @(k) EulerLogEq(h, x0, 12, r, k, t_obs);

%% Ground truth and Data
k_true = 12 * exp( - cos( pi/6 * ts) / 4 );
y_true = Gh(k_true);

% Noise variance
rng(17);
sigma = 0.5;
y = y_true + sigma * randn(1,12);

%% Plotting
figure();
ts = 0:h:12;
plot(ts, k_true,'-k','LineWidth',2);
hold on;
grid on;
plot(ts,  EulerLogEq(h, x0, 12, r, k_true), '-b','LineWidth',2);
plot(t_obs, y, 'ro','LineWidth',2);
xlabel('t')
set(gca,'fontsize',14)
legend('k_{true}', 'x_{true}', 'data')

%% Prior for log(k)

nt = length(ts);
N = 1000; % number of KL terms
KL = zeros(nt, N);
Lambda = zeros(N,1);

% Evaluating KL at grid points and setting Lambda 
KL(:,1) = ones(nt,1)/12;
Lambda(1) = 1/(pi^2);
for i = 2:N,
    KL(:,i) = sin((i-1) * pi/12 * ts)/sqrt(6); 
    Lambda(i) = 1/(pi * (i-1))^2;
end

% Covariance matrix and Cholesky factor for sampling
C0 = diag(Lambda);
L0 = sqrt(C0);

% function mapping KLE coefficients u to k
u2k = @(u) exp(log(10) + KL * u);

% Plotting 100 prior realizations and truth
figure();
plot(ts, u2k(L0 * randn(N,100)),'LineWidth',2); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title('Prior samples')
set(gca,'Fontsize',14)

%% Potential for KLE coefficients
Pot = @(u) norm(y - Gh(u2k(u))).^2/2/sigma^2;

%% Negative Log Prior Density
P0 = @(u) norm(u./sqrt(Lambda)).^2/2;

%% pCN-MCMC
M = 2e4; %lenght of chain
Us = zeros(N, M); % states of chain

% Quantity of interest
w = [0.5, ones(1, length(ts)-2), 0.5] / h; % trapezoidal rule weights
f = @(u) w * u2k(u);
F = zeros(1,M);

% Stepsize
s = 0.17; % pCN, N = 1000
s = 0.06; % RW, N = 1000

A = 0;
tic;
u = Us(:,1);
uPot = Pot(u); % pCN
uPot = Pot(u) + P0(u); % RW
for i = 1:M,
    % pCN
    %v = sqrt(1-s^2) * u + s * L0 * randn(N,1);
    %vPot = Pot(v);
    
    % RW
    v = u + s * L0 * randn(N,1); %RW
    vPot = Pot(v) + P0(v);
    
    a = min(1, exp(uPot - vPot));
    A = A + a;
    
    if rand(1) <= a,
        u = v;
        uPot = vPot;        
    end
    
    Us(:,i) = u;
    F(i) = f(u);
end
toc;
A/M

%% Plotting posterior

% 100 realizations and truth
figure();
plot(ts, u2k( Us(:,1:1e2:end)),'LineWidth',2); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title('Posterior draws')
set(gca,'Fontsize',14)

% Mean and 95 credible intervals
Ks = u2k(Us);
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
legend('Mean','5%-quantil','95%-quantil','truth','Location','NorthEast')
title('Posterior Mean and credible region')
set(gca,'Fontsize',14)

%

% Plot Path of Markov chain
figure();
plot(M-1000:M,F(M-1000:M),'o-','LineWidth',2);
grid on;
xlabel('Iteration')
set(gca,'FontSize',14)
title('Path of f(u)')


