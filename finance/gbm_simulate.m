% Sun 12 Jan 11:43:20 +08 2020
%
% S = gbm_simulate(t,r,sigma,S0,n)
%
% simulate geometric brownian motion
% 
% input :
% t : vector of times
% r : risk free interest rate
% sigma1 : volatility for unit time step
% S0 : start price
% n  : number of parallel simulations
function S = gbm_simulate(t,r,sigma1,S0,n)
	R  = randn(length(t)-1,n);
	dt = diff(t);
	
	W  = cumsum(dt.*R);
	S  = S0*[ones(1,n);
                 exp( r*(t(2:end)-t(1)) + sqrt(dt).*sigma1.*W)];
	% for skewed : 
%	sqrt(1-delta^2)*W1 + delta*abs(W2) 
end

