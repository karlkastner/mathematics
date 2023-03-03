% 2019-02-19 20:16:11.163129683 +0100
% Karl Kastner, Berlin
%
% cdf of the skewed log-normal distribution
function f = logskewcdf(x,mu,s,a)
	%f = normpdf((log(x)-mu)/s) - 2*normcdf((log(x)-mu)/s,a);
	f = zeros(size(x));
	for idx=1:prod(size((x)))
		f(idx) = quad(@(x_) logskewpdf(x_,mu,s,a),0,x(idx));
	end
end
