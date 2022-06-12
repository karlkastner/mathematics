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
	S  = spectral_density_lowpass(fx,rho,order,dx);
	y  = ifft(S.^(order/2).*fft(x));

% note : this is solely for test purposes, directly multiply the ft once
%        with the spectrum of the filter (!)
%	L = n-1;
%	D2 = fourier_derivative_matrix_1d(n,L,2);
%	rho = rho/(1 - 2*rho + rho^2);
	% (I - rho*D) \ y
%	%y = pcg(@fun, flat(x),[],100);
%	y = pcg(@fun, x,[],1000);
%
%	function y = fun(x)
%		f   = fft(x);
%		D2f = D2*f;
%		D2x = ifft(D2f);
%		y   = x - rho*D2x;
%	end
end

