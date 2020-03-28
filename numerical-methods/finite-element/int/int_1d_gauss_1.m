% Tue Feb 28 22:48:27 MSK 2012
% Karl Kästner, Berlin

% w : weights
%	2/(1-xi^2)(P'_n(xi))^2
% b : baricentric coordinates
%     ith-root of legendre polynomial of order n
% second order, midpoint rule
function [w, b, flag] = int_1d_gauss_1()
	% weights of integration points
	w = 1.0;
	% baricentric coordinates of integration points
	b = [0.5 0.5];
	% mark scheme to be diagonal
	flag = 0;
end

