% Sat 11 Jun 16:43:15 CEST 2022
function Ry = acf_brownian_phase_across(ly,sy)
	Ry = exp(-0.5*abs(ly).*sy.^2);
end

