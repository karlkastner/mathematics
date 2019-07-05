% Thu 21 Jul 11:07:40 CEST 2016
% Karl Kastner, Berlin
%
%% zeros of continuous fourier series series
%%
%%	f  = a_0 + sum_j=^n a_i cos(j x) + b_i sin(j x)
%
function [zF, zB, B] = zeros_fourier(a,b)
	n = length(a)-1;
	if (0 == n)
		% if the function is a constant, then there are no zeros
		zF = [];
		zB = [];
		B  = [];
	end

	if (isnumeric(a))
		h = zeros(2*n+1);
		B = zeros(2*n);
	end
	for k=0:n-1
		h(k+1) = a(n-k+1) + 1i*b(n-k+1);
	end
	h(n+1) = 2*a(1);
	for k=n+1:2*n
		h(k+1) = a(k-n+1) - 1i*b(k-n+1);
	end
	for k=0:2*n-1
		if (k>0)
			B(k+1,k) = 1;
		end
		B(k+1,2*n) = -h(k+1)/(a(n+1) - 1i*b(n+1));
	end
	% roots of B
	zB = eig(B);

	% roots of the fourier series
	zF = -1i*log(zB);
end

