% 2025-02-27 19:12:56.242462673 +0100
% 
%% determine correlation length of the process
%%	z = exp(lz)
%%	z  ~ lognormal(mu,sd)
%%	lz ~ normal(lmu,lsd)
%%	R_lz(x,y) = exp(-sqrt(x^2 + y^2)/lz)
%% correlation lengths
%%	Rz(t)   = exp(-1)
%%	Rzl(tl) = exp(-1)
%	
%
function t = geometric_ou_correlation_length_log2lin(ls,lt)
	t = -lt*log((log(exp(ls^2) + exp(1) - 1) - 1)/ls^2);
end

