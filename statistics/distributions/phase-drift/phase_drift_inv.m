function fx = brownian_phase_inv(cdf,f0,sx);
	a  = tan(pi*cdf);
	% TODO, this is not correct for for 
	error('account for sign switch')
%	fx = (f0*sqrt(a.^2*pi^2*sx^4 + a.^2 + pi^2*sx^4) - pi*f0*sx^2)./a;
	
	%a  = tan(pi*cdf-pi/2);
%	fx(:,2) = -(f0*sqrt(a.^2*pi^2*sx^4 + a.^2 + pi^2*sx^4) - pi*f0*sx^2)./a;
	fx = (sqrt(f0^2 - exp(cdf*pi*1i)) + f0*pi*sx^2*1i);
end

