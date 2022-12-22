
% note : this is solely for test purposes, directly multiply the ft once
function test_lowpass1d_fft(x,rho,order)
        with the spectrum of the filter (!)
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

