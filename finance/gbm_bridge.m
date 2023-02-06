% Mon 23 Jan 14:26:55 CET 2023
%
% S = gbm_bridge(t,r,sigma,S0,n)
%
% simulate geometric brownian motion
% 
% input :
% t : vector of times
% r : risk free interest rate
% sigma1 : volatility for unit time step
% S0 : start price
% n  : number of parallel simulations
function S = gbm_bridge(t,r,sigma1,S0,n)
	nt = length(t);
	e  = [0; randn(nt,n)];
	dt = diff(t);
	W  = cumsum(dt.*(e-0.5*e));

	% transform motion to bridge
	W  = W(1:end-1,:) - cvec((t-t(1))).*(W(end,:)-W(1,:));

	S  = S0*exp( r*(t(2:end)-t(1)) + sqrt(dt).*sigma1.*W);
	% for skewed : 
%	sqrt(1-delta^2)*W1 + delta*abs(W2) 
end

