% 2014-07-01 19:40:33 +0700
% Karl Kastner, Berlin
%
%% moments estimated from quantiles
%
% TODO, scale kr and sk similar to sigma, such that kr_q equals qr for sum of two normals
% for symmetric distribution input is expected to be |X-mu|
function [mu sigma sk kr] = qmoments(x,symmetric)
	if (nargin() < 2 || isempty(symmetric) || ~symmetric)
	% assymetric distribution
	% octiles of the distribution
	oct = quantile(x,(1:7)/8);
	% location
	% m1 = q2 = o4;
	mu = oct(4);
	% norminv(q50)

	% dispersion
	if (nargout() > 1)
		% q3 - q1
		sigma = (oct(6)-oct(2))/(2*0.6745);
		% (norminv(q75)-norminv(q25))/(2*0.6745)
	end
	% skewness
	if (nargout() > 2)
		% (bowleys 1920)
		% q1 - 2*q2 + q3
		sk = (oct(2) - 2*oct(4) + oct(6))/(oct(6)-oct(2));
	end
	% kurtosis
	if (nargout() > 3)
		% Moors (1988)
		kr = ((oct(7) - oct(5)) + (oct(3) - oct(1))) / (oct(6) - oct(2)) - 1.2331;
	end
	
	else
	% compute moments for symmetric distributions
	q = quantile(x,[0.25 0.5 0.75]);
	mu = 0;
	if (nargout()>1)
		sigma = q(2)/0.6745;
	end
	if (nargout()>2)
		sk = 0;
	end
	if (nargout()>3)
		kr = (q(3) - q(1))/q(2) - 1.2331;
	end
	
	end
end % function qmoments

