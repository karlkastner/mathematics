% 2016-02-12 19:06:55.520166486 +0100
% Karl Kastner, Berlin
%
%% coefficients of the ar2 process determined from the two leading correlations
%% of the acf [1,r1,r2,...]
% r1, r2: leading coefficients of the acf (sample estimates of rho1 and rho2)
% a1, a2: coefficients of the AR2-process Y_i - a1 Y_i-1 - a2 Y_i-2 = 0
function [a1 a2] = ar2conv(r1,r2)
	if (1 == nargin)
		r2 = r1(2);
		r1 = r1(1);
	end

	a1 = (r1*r2 - r1)/(1-r1^2);
	a2 = (r1^2  - r2)/(1-r1^2);

	if (nargout()<2)
		a1 = [a1 a2];
	end

