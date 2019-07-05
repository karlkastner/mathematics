% Wed 19 Jul 15:44:36 CEST 2017
% Karl Kastner, Berlin
%
%% variance of sub sample starting at the end of the series
%%
%% from the finite length first order autocorrelated process
%%
%% s2 = 1/m^2 sum_i^m sum_j^m rho^-|i-j|
%
function s2 = ar1_var_range2(sigma,rho,n,m)
%	ij    = bsxfun(@minus,(1:m)',(1:m));
%	rhoij = rho.^abs(ij);
%	sum_rhoij = sum(rhoij(:));
	sum_rhoij = m + 2*rho/(1-rho)*(m - (1-rho.^m)/(1-rho));
	s2 = sigma^2*(1-1./m.^2.*sum_rhoij);
end

