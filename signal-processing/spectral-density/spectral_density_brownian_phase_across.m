% Sat 11 Jun 15:06:15 CEST 2022
%
% function S = spectral_density_brownian_phase_across(fy,sy)
%
function S = spectral_density_brownian_phase_across(fy,sy)
	if (issym(fy) || issym(sy))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	ky = 2*pi_*fy;
	S = 4*sy.^2./(sy.^4+ky.*ky);
end

