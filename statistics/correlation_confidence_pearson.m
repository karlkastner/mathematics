% Thu 25 Aug 09:29:58 CEST 2016
% Karl Kastner, Berlin
%
%% confience intervals of the correlation coefficient
%% c.f. Fischer 1921
%
function rho_lu = correlation_confidence_pearson(rho,n,c)
	% transform to normal distributed value
	F  = 0.5*log((1+rho)/(1-rho));
	% standard error of the transformed value
	se = 1/sqrt(n-3);
	% the confidence intervals
	z  = 1-normcdf(c);
	% upper and lower bound for F
	F_lu = F + se*[-z,z];
%	rho_ = tanh(F_);
	% upper and lower bound for rho
	rho_lu = (exp(2*F_lu)-1)./(exp(2*F_lu)+1);
end
	
	
