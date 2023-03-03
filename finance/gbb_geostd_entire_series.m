function gsd = gbb_geostd_entire_series(t,S0,T,r,sigma)
	% varbar = s^2 T/6
	gsd = exp(sigma*sqrt(T/6));
end
