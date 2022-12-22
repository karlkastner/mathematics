% 2022-03-30 12:46:31.571547062 +0200
% TODO simplify to univariate optimization :
% xm = exp(mu - s^2)
% ym = 1./(xm*s*sqrt(2pi))*exp(- (ln(xm) - mu)^2/(2s^2))
% ym*xm = 1/(s*sqrt(2pi))*exp(-s^2/2)

function [lmu,lsd] = logn_mode2param(xm,ym)
	l = lsqnonlin(@res,[log(xm),xm]);
	lmu = l(1);
	lsd = l(2);

function res = res(l)
	[xm_,ym_] = logn_mode(l(1),l(2));
	res = [xm-xm_; ym-ym_];
end

end

