% 2016-03-03 13:18:19.590258695 +0100
% Karl Kastner, Berlin
%
%% log normal distribution called by modes rather than parameters
function f = lognpdf_(x,mu,sd)
	% convert parameters
	[lmu lsd] = logn_mode2param(mu,sd);
	%[lmu lsd] = logn_param2mode(mu,sd);
	
	% evaluate
	f = lognpdf(x,lmu,lsd);
end

