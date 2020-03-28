% 2016-03-03 13:52:20.475878559 +0100
% Karl Kastner, Berlin
%
% polynomial fit to pdf
%
function [pdf, x] = polypdf(x,x0,order)
	x0    = sort(cvec(x0));
	n     = length(x0);
	cdf0  = (1:n)'/(n+1);

	A    = vander_1d(x0,order);
	ccdf = A\cdf0;

	if (isempty(x))
		q = quantile(x0,0.01,0.99);
		x    = linspace(q(1),q(2),100)';
	end
	A    = vander_1d(x,order);
	cpdf = cvec(polyder(flipud(ccdf)));
	pdf  = A*ccdf;
end

