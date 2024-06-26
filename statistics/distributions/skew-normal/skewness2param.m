% 2015-08-03 14:09:44.313605435 +0200
% Karl Kastner, Berlin
%
% convert skewness to skew-parameter of the skew-normal-distribution
% azzalini
function a = skewness2param(sk)
	lim = skewnormal_maxsk();
	if (abs(sk) > lim)
		warning(['skewness has to be between +/-',num2str(lim)]);
	end
	sk23  = abs(sk).^(2/3);
	%delta = sqrt(sk23/(sk23*2/pi + (0.5*(4-pi))^(2/3)*2/pi));
	delta = sqrt(sk23./(2/pi*(sk23 + (2-pi/2)^(2/3))));
	a     = sign(sk).*delta./sqrt(1-delta.^2);
end

