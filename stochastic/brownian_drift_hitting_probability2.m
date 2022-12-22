% Mon  4 Jul 16:10:49 CEST 2022
% c.f. abundo
% p tau < x
function P = p_passage(t,a,b,x)
	P = 1 - normcdf((a-x)./sqrt(t)+b.*sqrt(t)) ...
	      + exp(-2*b*(a-x)).*normcdf(-(a-x)./sqrt(t)+b.*sqrt(t));
end

