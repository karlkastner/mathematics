% Fri 12 Oct 12:42:33 CEST 2018
% Karl Kastner, Berlin
%
%% derivative by fourier transform
%% exponential convergence for periodic functions
%% results in spurious oscillations for aperiodic functions
%%
%% input:
%% x : data, sampled in equal intervals 
%% k : order of the derivative
%%
%% dx : kth-derivative of x
%
%
% TODO what about interval length? 1/dx^k
% TODO flag for not computing the fft
function dx = fft_derivative(x,k)
	if (nargin()<2)
		k = 1;
	end
	x = cvec(x);
	x = fft(x);
	n = length(x);
	if (1 == mod(n,2))
		m = (n-1)/2;
		id = (1:m)';
		id_= (m:-1:1)';
		dx = (2i*pi/n).^k*[0;
                                     id.*x(2:m+1);
                                    -id_.*x(m+2:2*m+1)];
	else
		m  = n/2-1;
		id = (1:m)';
		id_= (m:-1:1)';
		dx = (2i*pi/n).^k*[0;
                                   id.*x(2:m+1);
			           0;
                                   -id_.*x(m+3:2*m+2)]/n;
	end
	dx = ifft(dx);
end

