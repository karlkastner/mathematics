% 2021-09-23 10:17:05.746154118 +0200
% flat spectral density of a random vector with iid elements
function S = spectral_density_iid(L,n)
	S = 2*L/n*ones(n,1);
end

