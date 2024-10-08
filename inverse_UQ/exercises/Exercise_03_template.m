% Template for Exercise 03: Markov chain Monte Carlo simulation

clear all; clc; close all;

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
...
    
% Covariance matrix and Cholesky factor for sampling
...
    
% function mapping KLE coefficients u to k
u2k = @(u) ...
    
% Plotting 100 prior realizations and truth
figure();
plot(ts, ... ,'LineWidth',2); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title('Prior samples')
set(gca,'Fontsize',14)

%% Potential for KLE coefficients
Pot = @(u) ...

%% Negative Log Prior Density
P0 = @(u) ...

%% MCMC
M = ... ; %lenght of chain
Us = zeros(N, M); % states of chain

% Quantity of interest
w = ... % trapezoidal rule weights
f = @(u) ... 
F = zeros(1,M); % values of f at states of chain

% Stepsize
s = ...

% average acceptance rate
A = 0;
tic;
u = Us(:,1);
uPot = ... % potential at current state
for i = 1:M,
    % proposal and its potential
    ...
    
    
    % Computing acceptance probability
    ...

    % Accept reject step
    ...
    
    % Storing current state
    Us(:,i) = ...
    F(i) = ...
end
toc;

% average acceptance rate
A/M

%% Plotting posterior

% 100 realizations and truth
figure();
plot(ts, u2k( Us(:,1:1e2:end) ),'LineWidth',2); hold on
plot(ts, k_true,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
title('Posterior draws')
set(gca,'Fontsize',14)

% Mean and 95 credible intervals
...

figure();
plot(ts, ...,'-b','LineWidth',2); hold on
plot(ts, ...,'--r','LineWidth',2);
plot(ts, ...,'--r','LineWidth',2);
plot(ts, ...,'-k','LineWidth',2);
xlabel('t')
ylabel('k(t)')
grid on;
legend('Mean','5%-quantil','95%-quantil','truth','Location','NorthEast')
title('Posterior Mean and credible region')
set(gca,'Fontsize',14)

% Compute asymptotic confidence interval for posterior mean of f
[acf] = autocorr(F,'NumLags',...)
...

