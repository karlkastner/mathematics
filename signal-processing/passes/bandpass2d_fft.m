% Sun 11 Jul 22:56:09 CEST 2021
% note : this function is for testing purposes only,
%        directly multiply the ft of the signal with the ft of the filter
%        to obtain the filtered signal in a single step
%
% function y = bandpass2d_fft(x,rho,a,order)
function y = bandpass2d_fft(x,rho,a,order)
	if (nargin()<4)
		order = 1;
	end

	if (length(rho) == 1)
		rho = rho*[1,1];
	end

	n = size(x);
	L = n-1;
	[Dx, Dy, D2x, Dxy, D2y] = fourier_derivative_matrix_2d(n,L);

	% renormalize
	rho = rho./(1 - 2*rho + rho.*rho);

	% derivative matrix
	if (nargin()<2 || isempty(a))
		% no rotation, coordinate axis-parallel smoothing
		D2 = rho(1)*D2x+rho(2)*D2y;
	else
		c = cos(a);
		s = sin(a);

		% derivative along diangonal with angle a
		D2s = c*c*D2x + 2*s*c*Dxy + s*s*D2y;
		% derivative along orthogonal direction
		D2n = s*s*D2x - 2*s*c*Dxy + s*s*D2y;
		% combined
		D2 = rho(1)*D2s + rho(2)*D2n;
	end
	
	% y = (I-r*D2)^-2*D2*x
	% note that n-iterations with D2 is faster
	% than a single step with D2^n, as D2^n is less well conditioned
	for idx=1:order

	% high pass, derive with D2
	% y = D2*x
	x = flat(x);
%	f   = flat(fft2(x));
%	D2f = D2*f;
%	D2f = reshape(D2f,n);
%	x   = flat(ifft2(D2f));
	x = fun(x);	

	% low pass
	% y = (I - rho*D)*x
	x = pcg(@fun, flat(x),[],sum(n));
	x = reshape(x,n);
	end

	y = x;

	function y = fun(x)
		x_  = reshape(x,n);
		f   = flat(fft2(x_));
		D2f = D2*f;
		D2f = reshape(D2f,n);
		D2x = flat(ifft2(D2f));
		y   = x - D2x;
	end % fun

%	function y = fun4(x)
%		x_  = reshape(x,n);
%		f   = flat(fft2(x_));
%		D2f = (I -2*D2 + D4)*f;
%		D2f = reshape(D2f,n);
%		D2x = flat(ifft2(D2f));
%		y   = D2x;
%	end
end % bandpass2d_fft

