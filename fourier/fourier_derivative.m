% Tue  7 Aug 14:03:50 CEST 2018
% Karl Kastner, Berlin
%
%% coefficients of the derivative of a fourier series
%% not of discrete fourier transform (fft) 
%
% TODO move all fourier-series functions into the class fourier-series
function [d] = fourier_derivative(T, c)
	d = zeros(size(c));
	omega = 2*pi./cvec(T);
	d(2:2:end-1,:) =  bsxfun(@times,omega(1:end),-c(3:2:end,:));
	d(3:2:end,:)   =  bsxfun(@times,omega(1:end),+c(2:2:end-1,:));
%	t   = t-t0(1);
%	F   = fourier_matrix(T, t);
%	val = F*c;

%	A = zeros(nt,2*no+1);
%	A(:,1) = 1;
%	for idx=1:no
%		A(:,2*idx)   = sin(2*pi*t/T(idx));
%		A(:,2*idx+1) = cos(2*pi*t/T(idx));
%	end
%	val = A*c; %cvec(c);
	% TODO prediction error
end

