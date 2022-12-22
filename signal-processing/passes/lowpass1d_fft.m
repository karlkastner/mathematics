% Sun 11 Jul 21:39:38 CEST 2021
%
% function y = lowpass1d_fft(x,rho,order)
function y = lowpass1d_fft(x,rho,order)
	if (nargin()<3)
		order = 1;
	end
	% order_fft = sqrt(order_implicit)
	order = sqrt(order);
	if (isvector(x))
		x = cvec(x);
	end
	n  = size(x,1);
	fx = fourier_axis(1,n);
	dx = 1/n;
	S  = spectral_density_lowpass_discrete(fx,rho,order,dx);
	y  = ifft(S.^(order/2).*fft(x));

end

