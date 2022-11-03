% Tue 21 Sep 13:24:55 CEST 2021
% m : number of frequencies bins the max Shat was chosen from
function p = periodogram_p_value(Shat,S,m)
	dof = 2;
	if (nargin()<3)
		m = 1;
	end
	p1 = 1-chi2cdf(dof*Shat/S,dof);
	% correct for multiple testing
	% c = Shat/S ~ 1/2*chi2inv((1-p)^(1/m),2)
	% p = 1-chi2cdf(2*Shat/S)^m
	p  = 1-(1-p1)^m;
end

