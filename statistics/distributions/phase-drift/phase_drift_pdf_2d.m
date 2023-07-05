function [S,Sx,Sy] = spectral_density_brownian_phase_2d(fx,fy,f0,sx,sy,normalize)
	if (nargin()<6)
		normalize = true;
	end
	if (isvector(fx) && isvector(fy))
		fx = cvec(fx);
		fy = rvec(fy);
	end

	% perpendicular to bands
	Sx = spectral_density_brownian_phase(fx,f0,sx,false);
	% parallel to bands
	Sy = spectral_density_brownian_phase_across(fy,sy);
	% 2D density
	S   = Sx*Sy;
	if (normalize)
		dfx = fx(2)-fx(1);
		dfy = fy(2)-fy(1);
		% normalize
		S = 2*S/(sum(sum(S))*dfx*dfy);
	end
end

