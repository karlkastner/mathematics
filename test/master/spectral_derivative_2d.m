% Sat Aug  4 13:58:28 MSK 2012
% Karl Kästner, Berlin

% spectral derivative matrices
% u must be shaped into an 2d array
% expects even number of input samples in both directions
% expects lowest frequencies on boundary (as returned by fft)
function [Dx Dy] = df_2d(n, p)
%function [du_dx du_dy] = df_2d(n, p)
	n = n/2;
	% rows (x)
%	Dx = [2*n(1):-1:1]';
%	Dx = [+n(1) -(n(1):-1:1) +(2:+1:n(1))]';
%	Dx = [+n(1) -(n(1)-1:-1:0) +(1:+1:n(1)-1)]';
%	Dy = -[+n(2) -(n(2):-1:1) +(2:+1:n(2))];
	Dx = -[+(0:+1:n(1)) -(n(1)-1:-1:1)]';
	Dx = ((-1i*2*pi*Dx).^p)*ones(1,2*n(2));
	% columns (y)
	Dy = -[+(0:+1:n(2)) -(n(2)-1:-1:1)];
	Dy = ones(2*n(1),1)*((1i*2*pi*Dy).^p);
end % df_2d

