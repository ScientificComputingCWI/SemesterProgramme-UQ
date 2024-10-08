function x = EulerLogEq(h, x0, T, r,  k, t_obs)

% This function implements the Euler scheme for the logistic equation with 
% time-dependent capacity
% dx/dt = r x(t) (k(t) - x(t))
% with initial condition
% x(0) = x0
%
% Input:
% h ... stepsize of discretization (scalar)
% x0 ... initial state (column vector)
% T ... end point in time (scalar)
% r ... reproduction rate
% k ... capacity (column vector of correct length)
% t_obs ... times of observation (optional and multiple of h)
%
% Output:
% x ... matrix consisting of states x at time k*h, k in N, or at
%       observation times

% vector of times and solution
t = 0:h:T;
nt = length(t);
x = zeros(1, nt);

% Initialization
x(1) = x0;

for i = 2:nt,
    x(i) = x(i-1) + h * r * x(i-1) * (k(i-1) - x(i-1)); 
end
    
if nargin == 6,
    ind = 1 + round(t_obs/h);
    x = x(ind);
end