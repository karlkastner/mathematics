% 2015-10-15 16:43:10.980502041 +0200
% Karl Kastner, Berlin
%
%% round to n digits
%
function X = roundn(X,digits)
	fdx = (0==X);
	E = ceil(log10(abs(X)));
	X = X.*10.^-E;
	X = round(X,digits);
	X = X.*10.^E;
	X(fdx) = 0;
end

