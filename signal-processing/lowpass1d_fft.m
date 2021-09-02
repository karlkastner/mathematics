% Sun 11 Jul 21:39:38 CEST 2021
% note : this is solely for test purposes, directly multiply the ft once
%        with the spectrum of the filter (!)
function y = lp1d_fft(x,rho)
	n = length(x);
	L = n-1;
	D2 = fourier_derivative_matrix_1d(n,L,2);

	rho = rho/(1 - 2*rho + rho^2);
	% (I - rho*D) \ y


	%y = pcg(@fun, flat(x),[],100);
	y = pcg(@fun, x,[],1000);

	function y = fun(x)
		f   = fft(x);
		D2f = D2*f;
		D2x = ifft(D2f);
		y   = x - rho*D2x;
	end
end

