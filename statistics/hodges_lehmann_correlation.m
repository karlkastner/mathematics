% Thu 27 Sep 16:47:18 CEST 2018
%
%% hodges_lehmann correlatoon coefficient
%% c.f. Shamos 1976
%% c.f. Bickel and Lehmann 1976
%% c.f. rousseeuw 1993
%% c.f. Shevlyakov 2011
function rho = hodges_lehmann_correlation(x,y)
	if (length(x) ~= length(y))
		error('here');
	end
	
	% normalise
	x = x/hodges_lehmann_dispersion(x);
	y = y/hodges_lehmann_dispersion(y);
	
	% identitiy of correlation coefficient
	u = x+y;
	v = x-y;
	s2u = hodges_lehmann_dispersion(u)^2;
	s2v = hodges_lehmann_dispersion(v)^2;
	rho = (s2u-s2v)/(s2u+s2v);

end

