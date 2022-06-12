% Tue  7 Aug 14:03:50 CEST 2018
% Fri 12 Oct 12:42:33 CEST 2018
% Karl Kastner, Berlin
%
%% derivative via fourier transform
%% exponential convergence for periodic functions
%% results in spurious oscillations for aperiodic functions
%%
%% input:
%% x : data, sampled in equal intervals 
%% k : order of the derivative
%%
%% dx : kth-derivative of x
%%
%% note : 1) the derivative converges with spectral accuracy, i.e. is
%%           exact up to rounding condition for L sufficiently large
%%	     and x being periodic
%%	  2) the derivative converges with order p, when x has only
%%	     p-continous derivatives, including discontinuous derivatives
%%	     over the boundary
%%	  3) discontinuous derivatives result in gibbs phenomenon
function dx = fft_derivative2(x,L,k)
	if (nargin()<2||isempty(L))
		L = 1;
	end
	if (nargin()<3)
		k = 1;
	end
	if (isvector(x))
		x = cvec(x);
	end
	x = fft(x);
	n = size(x,1);
%	if (1 == mod(n,2))
%		m = (n-1)/2;
%		%id = (1:m)';
%		%id_= (m:-1:1)';
%		%dx = (2i*pi/n).^k*[0;
%                %                     id.*x(2:m+1);
%                %                    -id_.*x(m+2:2*m+1)];
%		%o = (2i*pi*(1:m)'/n);
%		%o_= (-2i*pi*(m:-1:1)'/n);
%		%dx = [0;
%                %                     o.^k.*x(2:m+1);
%                %                    o_.^k.*x(m+2:2*m+1)];
%		%f = [ (0:m)';
%		%     -(m:-1:1)']/n;
%		f = fourier_axis(n,n)';
%		dx = (2i*pi*f).^k.*x;
%	else
%		m  = n/2-1;
%		id = (1:m)';
%		id_= (m:-1:1)';
%		dx = (2i*pi/n).^k*[0;
%                                   id.*x(2:m+1);
%			           -(m+1).*x(m+2);
%                                   -id_.*x(m+3:2*m+2)];
%	end
	fx = cvec(fourier_axis(L,n));
	dx = (2i*pi*fx).^k.*x;
	dx = ifft(dx);
end

