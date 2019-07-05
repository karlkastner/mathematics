% 2015-01-07 13:54:51.401372485 +0100
% Karl Kastner, Berlin
%
%% transformation matrix for a continuous fourier series
%% (not for the discrete dft/fft)
%
function A = fourier_matrix2(n,m,X,L)
	if (isempty(X))
		X = 2*pi*(0:n-1)'/n;
	else
		if (nargin()<4)
			L = [min(X) max(X)];
		end
		X = 2*pi*(X-L(1))/(L(2)-L(1));
		n = length(X);
	end
	A = zeros(n,2*m+1);
	A(:,1) = 1;
	for idx=1:m
		A(:,2*idx)   = sin(idx*X);
		A(:,2*idx+1) = cos(idx*X);
	end
end

