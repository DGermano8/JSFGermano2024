function [Z,t] = TauLeapingMethod(x0, rates, stoich, times, options)
%% The Tau Leaping Method

% initialise

X0 = x0;
nu = stoich.nu;

tFinal = times(end);

try tau = options.dt; catch tau = 10^(-3);
end

Nt = floor(tFinal/tau);
Z = zeros(length(X0),Nt+1);
t = zeros(1,Nt+1);
Z(:,1) = X0;

for i=1:Nt
    % compute propensities
    Props = rates(Z(:,i),t(i));
    Props(Props < 0 ) = 0;
    % generate poisson variates
    Y = poissrnd(Props*tau);
    % update copy numbers
    Z(:,i+1) = Z(:,i) + (nu') * Y;
    t(i+1) = t(i) + tau;
end
