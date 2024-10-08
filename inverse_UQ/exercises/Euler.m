function x = Euler(h, x0, T, f, t_obs)
% This function implements the Euler scheme for an ODE
% dx/dt = f(t,x)
% with initial condition
% x(0) = x0
%
% Input:
% h ... stepsize of discretization (scalar)
% x0 ... initial state (column vector)
% T ... end point in time (scalar)
% f ... RHS of ODE (function handle for arguements t and x)
% t_obs ... times of observation (optional and multiple of h)
%
% Output:
% x ... matrix consisting of states x at time k*h, k in N, or at
%       observation times

% vector of times and solution
t = 0:h:T;
nt = length(t);
x = zeros(length(x0), nt);

% Initialization
x(:,1) = x0;

for i = 2:nt,
    x(:, i) = x(:, i-1) + h * f(t(i-1), x(:,i-1)); 
end
    
if nargin == 5,
    ind = 1 + round(t_obs/h);
    x = x(:,ind);
end
