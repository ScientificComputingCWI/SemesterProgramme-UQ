%
%********************************************************************
%***  Elisabeth Ullmann                                        ******
%***  Autumn School on Uncertainty Quantification @ CWI        ******
%********************************************************************
%
% Cross-entropy method with multivariate Gaussian density to estimate 
%
% P_f = Prob (min(X)>a) = int_{-min(x)+a<0} f(x)dx,
%
% where f(x) ~ N(0,Id) and a>0
%
clear all; 

a = 2;
Pf_e = normcdf(-a)^2;    % Exact probability of failure

d=2;                     % Input dimension
N = 1e3;                 % Number of samples
rho = 0.1;               % Quantile related to failure level

mu0=zeros(d,1);
Sigma0=[1 0; 0 1];
mu_ref=mu0;
Sigma_ref=Sigma0;
 
MU_ref = [mu_ref];       % Matrix to store auxiliary means in CE method

G = @(x)(-min(x)+a);
s = 1e-6;               % Scaling parameter to regularize Sigma_ref

rng('default')           % For reproducibility
iter = 0;
f_level = 1;

while f_level>0

    x = mvnrnd(mu_ref',Sigma_ref,N)';
    y = G(x);
    f_level = quantile(y,rho);
    if f_level<0, break; end

    W = mvnpdf(x',mu0',Sigma0)./mvnpdf(x',mu_ref',Sigma_ref);
    H = W';
    ind = find(y<f_level);
    denom = sum(H(ind));

    % Update mean of Gaussian IS density
    mu_ref = sum(repmat(H(ind),d,1).*x(:,ind),2)/denom;

    % Update variance-covariance matrix of Gaussian IS density
    Sigma_ref = zeros(size(Sigma0));
    for k=1:length(ind)
        v_ref = x(:,ind(k))-mu_ref;
        Sigma_ref = H(ind(k))*v_ref*v_ref';
    end
    s = max(s,min(W));   % Adjust scaling parameter
    Sigma_ref = Sigma_ref/denom+s*eye(d);

    iter = iter+1;
    N_eff(iter) = sum(W)^2/sum(W.^2);
    MU_ref = [MU_ref mu_ref];
end
% Perform IS with final Gaussian density
W = mvnpdf(x',mu0',Sigma0)./mvnpdf(x',mu_ref',Sigma_ref);
N_eff(iter+1) = sum(W)^2/sum(W.^2);
H = W;
ind = find(y<0.0);
% CE-IS estimate of Pf
Pf_est_ce = sum(H(ind))/N;

% Standard Monte Carlo estimate of Pf
x = mvnrnd(mu0',Sigma0,N)';
y = G(x);
ind = find(y<0.0);
Pf_est_mc = length(ind)/N;

% Print results
fprintf('CE-estimate of Pf: %1.6e\n',Pf_est_ce)
fprintf('MC estimate of Pf: %1.6e\n',Pf_est_mc)
fprintf('Exact Pf: %1.6e\n',Pf_e)
fprintf('Number of CE iterations: %.0f\n',iter)
fprintf('Number of samples per level: %.0f\n',N)
% 
figure
plot(MU_ref(1,:), MU_ref(2,:),'ro','LineWidth',2)
grid on
xlabel('$x_1$',Interpreter='latex')
ylabel('$x_2$',Interpreter='latex')
title(sprintf('Mean of Gaussian CE approximation for a=%.0f',a))
