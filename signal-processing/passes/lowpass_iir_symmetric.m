% 2018-02-02 13:35:50.230627639 +0100
%% two-sided iir low pass filter (for symmetry)
function x = lp_iir_twosided(x,rho);
	if (isvector(x))
		x = cvec(x);
	end
	scale = (1-rho);
	if (0) % two-sided, fourth order
		x = filter(scale,[1,-rho],x,x(1));
		x = flipud(x);
		x = filter(scale,[1,-rho],x,x(1));
		x = flipud(x);
	else  % second order, one sided
		x = 0.5*(  filter(scale,[1,-rho],x,x(1)) ...
			 + flipud(filter(scale,[1,-rho],flipud(x),x(end))) ...
			 );
	end
end

