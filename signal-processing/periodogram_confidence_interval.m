% Fri 16 Jul 22:03:44 CEST 2021
% Karl Kästner, Berlin
%
%% confidence interval for periodogram values
%
% m : number of frequencies bins the max-value was choses from
function c  = periodogram_confidence_interval(p,S,m)
	dof = 2;
	p_  = (1-(1-p).^(1./m));
	c   = 1/2*chi2inv(p,dof).*S;
end

