% Do 4. Feb 10:14:20 CET 2016
% Karl Kastner, Berlin
%
%% scale factor for the standard deviation
%% of the asymtpotic distibution of sample quantiles
%% (for normal distribution)
%% see cadwell, 1952
function s = scale_sd_quant(p)
	f  = @normpdf;
	Fi = @norminv;
	q  = 1-p;
	c  = f(Fi(p));
	s  = sqrt(p.*q)./c;
end

