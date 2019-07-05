% Tue 19 Feb 20:36:56 CET 2019
%
% pdf of the skewed log distribution
function f = logskewpdf(x,mu,s,a)
	f = 2./(s*x).*normpdf((log(x)-mu)/s).*normcdf(a*(log(x)-mu)/s);
end

