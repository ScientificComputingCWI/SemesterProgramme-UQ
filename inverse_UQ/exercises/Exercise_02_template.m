% Template for Exercise 02: Approximation of Posterior

clear all; clc; close all;

%% (I) Forward map and its discretization
%
r = 0.25;
x0 = 0.1;

% Time of observation and their corresponding indices 
t_obs = [1,2,3];

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
Pot = @(k) ...
    
% Prior density
pi0_fun = @(k) 1/10;

% Computing posterior density at grid points ks
piy = zeros(length(ks),1);
...

% Computing normalizing constant Z via trapezoidal rule
% quadrature weights
w = ...
Z = ...
piy = piy/Z; % normalized posterior

figure(1);
plot(ks, piy,'-k'); 
title('Posteriors')
xlabel('k')
hold on;
grid on;
names{1} = 'true'
set(gca,'Fontsize',14)

%% Approximation of posterior

% Discretization parameters
hs = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];

% Approximate Potential
Poth = @(k,h) ...

dHell = zeros(1, length(hs)); % vector of Hellinger distances

% Loop over values for h
for j = 1:length(hs),
    
    ...
    
    % Plotting approximate posterior for certain h (where j is odd)
    if mod(j,2) == 1,
        figure(1);
        plot( ks, ... , '--'); 
            end

    % Computing Hellinger distance
    dHell(j) = ...
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
figure(2);
loglog(hs, dHell,'-ob')
xlabel('h')
ylabel('Helinger distance')
title('Error in Hellinger distance')
grid on
set(gca,'Fontsize',14)