% Fri Jan 30 16:57:03 CET 2015
% Karl Kastner, Berlin
%
%% transformation matrix for a continuous fourier series
%% (not for the discrete dft/fft)
%
% A = fourier_matrix(T,t)
function A = fourier_matrix_exp(T,t)
	% usually the first fourier coefficient has the pre-factor 1/2
	% here it is left with 1 so that it equals the mean (in case frequencies are orthogonal)
	no     = 1+2*length(T);
	A      = zeros(length(t),no);
	A(:,1) = 1;
	for idx=1:length(T)
		% was 2exp
		A(:,2*idx)   = exp(+2i*pi*t/T(idx));
		A(:,2*idx+1) = exp(-2i*pi*t/T(idx));
	end % for
end % fourier_matrix_exp

