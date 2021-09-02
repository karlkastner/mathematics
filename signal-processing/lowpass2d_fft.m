% Sun 11 Jul 21:39:38 CEST 2021
% note : this function is for testing purposes only,
%        directly multiply the ft of the signal with the ft of the filter
%        to obtain the filtered signal in a single step
function y = lowpass2d_fft(x,rho,a,order)
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
	if (nargin()<2)
		% no rotation, coordinate axis-parallel smoothing
		D2 = rho(1)*D2x+rho(2)*D2y;
	else
		c = cos(a);
		s = sin(a);

		% derivative along diangonal with angle a
		D2s = c*c*D2x + 2*s*c*Dxy + s*s*D2y;

		% derivative in orthogonal direction
		D2n = s*s*D2x - 2*s*c*Dxy + c*c*D2y;

		D2 = rho(1)*D2s + rho(2)*D2n;
	end

	y = flat(x);
	for idx=1:order
	end

		% y = (I - rho*D)*x
		y = pcg(@fun, y,[],sum(n));
		y = reshape(y,n);
	end

	function y = fun(x)
		x_  = reshape(x,n);
		f   = flat(fft2(x_));
		D2f = D2*f;
		D2f = reshape(D2f,n);
		D2x = flat(ifft2(D2f));
		y   = x - D2x;
	end % fun
end % lowpass2d_fft

