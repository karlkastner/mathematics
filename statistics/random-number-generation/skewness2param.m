% 2015-08-03 14:09:44.313605435 +0200
% Karl Kastner, Berlin
function a = skewness2param(sk)
	if (abs(sk) > 1-sqrt(eps))
		error('skewness has to be between -1 and 1');
	end
	sk23 = abs(sk).^(2/3);
	delta = sqrt(sk23/(sk23*2/pi + (0.5*(4-pi))^(2/3)*2/pi));
	a = sign(sk)*delta/sqrt(1-delta^2);
end

