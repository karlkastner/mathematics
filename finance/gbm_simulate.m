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
	if (nargin()<5)
		n = 1;
	end
	t  = cvec(t);
	nt = length(t);
	e  = randn(nt,n);
	dt = diff(t);

	% that is a stochastic integral, so the step is sqrt(dt)!
	W = [zeros(1,size(e,2)); 0.5*cumsum((e(1:end-1,:)+e(2:end,:)).*sqrt(dt))];

	S  = S0*exp( r*(t-t(1)) + sigma1.*W);
	% for skewed : 
%	sqrt(1-delta^2)*W1 + delta*abs(W2) 
end

