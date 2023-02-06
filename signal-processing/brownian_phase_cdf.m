% Thu 19 Jan 10:01:45 CET 2023
function cdf = brownian_phase_cdf(fx,f0,sx)
	if (~issym(fx))
		pi_ = pi;
	else
		pi_ = sym('pi');
	end

	cdf = atan2(2*f0*fx*pi_*sx^2,pi_^2*f0^2*sx^4 + f0^2 - fx.^2)/pi_;
end


