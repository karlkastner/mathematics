% 2016-05-09 18:44:11.290721154 +0200
% Karl Kastner, Berlin
%
%% floor to leading n-digits
%
function X = ceiln(X,d)
	% was round
	E = floor(log10(abs(X)));
	X = X.*10.^(-E+d-1);
	X = ceil(X);
	X = X.*10.^(E-d+1);
end

