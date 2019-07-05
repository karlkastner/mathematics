% Do 25. Feb 14:12:19 CET 2016
% Karl Kastner, Berlin
%
%% parameter covariance of the least squares regression
%%
%% fun : model function for predtiction
%% b   : sample values
%% f(p) = b
%% p   : parameter at point of evaluation (preferably optimum)
function [sp, C] = lsq_sparam(fun,b,p)
	% for OLS : H = 2*A'A, g = 2*A'b
	[H, g, f] = hessian(@(p) sum((fun(p) - b).^2)  ,cvec(p));
	% number of samples
	n = length(b);
	% number of parameter 
	np = length(n);
	% standard error
	serr2 = f/(n-np);
	% parameter covariance matrix
	C = serr2*inv(0.5*H);
	% standard error of parameters
	sp = sqrt(diag(C));
end

