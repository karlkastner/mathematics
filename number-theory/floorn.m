% 2016-05-09 18:44:11.290721154 +0200
% Karl Kastner, Berlin
%
%% floor to n-digits
%
function X = floorn(X,d)
	E = floor(log10(abs(X))) % was round
	X = X.*10.^(-E+d-1);
	X = floor(X);
	X = X.*10.^(E-d+1);
end

