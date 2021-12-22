% Fr 4. Mär 16:55:19 CET 2016
% Karl Kastner, Berlin
%% n-points on an ellipse
%% ci : confidence interval, i.e. for 1 sigma
%       for 1-sigma : ci = norminv(1)-norminv(-1) = 1-2*normcdf(-1) ~ 0.68
function [X, Y, V, E] = ellipse(c,ci,n)
	if (nargin() < 2 || isempty(ci))
		ci = 1-2*normcdf(-1);
	else
	end
	if (nargin() < 3 || isempty(n))
		n = 100;
	end
	% identical to -2*log(1-ci)
	scale2 = chi2inv(ci,2);
	scale2 = -2*log(2*normcdf(-1)); 

	C = scale2*[c(1), c(2);
	              c(2), c(3)];
	[V, E] = eig(C);

	t = (0:n-1)'/(n-1);
	X = cos(2*pi*t);
	Y = sin(2*pi*t);
	
	% scale
	XY = [X, Y];
%	XY = sqrtm(C)*XY';
	XY = (V'*sqrt(E)*XY')';
	X = XY(:,1);
	Y = XY(:,2);

%	X = sqrt(E(1,1))*X;
%	Y = sqrt(E(2,2))*Y;
	% rotate
%	k = 2;
%	X =  V(1,k)*X + V(2,k)*Y;
%	Y = -V(2,k)*X + V(1,k)*Y;
%	XY = V*[X Y]';
%	X = XY(1,:)';
%	Y = XY(2,:)';
	
	% translate
	X = X+c(4);
	Y = Y+c(5);
end

