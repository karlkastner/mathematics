% Sat 11 Jun 16:45:45 CEST 2022
function [R,Rx,Ry] = acf_brownian_phase_2d(lx,ly,fc,sx,sy)
	if (isvector(lx) && isvector(ly))
		lx = cvec(lx);
		ly = rvec(ly);
	end
	% perpendicular to bands
	Rx = autocorr_brownian_phase(lx,fc,sx);
	% parallel to bands
	Ry = autocorr_brownian_phase_across(ly,sy);

	R = Rx.*Ry;
end

