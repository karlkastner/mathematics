% Fri Oct 24 16:29:16 CEST 2014
% Karl Kastner, Berlin
%
%% fit a linear function by splitting the x-values at their median
%% 	(med(y_left) - med(y_right))/(med(x_left)-med(x_right)
%% this approach performs poorly compared to the theil-senn operator
%
% uncertainty in the prediction is
% s_p^2 = s^2_mu + (x-mu_x)^2*s^2_slope
%
% function [intercept slope var_e s_inter s_slope R2] = robustlinreg(X,Y,q)
function [c, var_e, s_inter, s_slope, R2] = robustlinreg(X,Y,q)
error('this seems to be broken, use medianslope')
	if (nargin() < 3 || isempty(q))
		q = 0.5;
	end
%	p = normcdf(-1);
	p = 0.25;
	w = -0.5*norminv(p);

	mx        = median(X);
	ldx       = X <= mx;
	rdx       = X >  mx;
	myl       = quantile(Y(ldx,:),q);
	mxl       = median(X(ldx));
	myr       = quantile(Y(rdx,:),q);
	mxr       = median(X(rdx));
	slope     = ( myr - myl ) / ( mxr - mxl );
	% both left and right lead to the same slope estimate,
	% however, the "average" is numrerically more stable (really? should depend on sign!) 
	% note : using intercept based on all values here does not yield the most efficient slope
	intercept  = 0.5*(myr+myl) - 0.5*(mxr+mxl)*slope;
	% residual
	res = (intercept + X*slope) - Y;
	n   = length(X);
	if (1)
		% error estimate based on moments
		% variance of the residuals (aka mse)
		var_e = (res'*res)/(n-2);
		var_y = var(Y);
		var_x = var(X);
	else
		% error estimate based on quantiles

		% variance of the residual
		qe   = quantile(res,[p 1-p]);
		var_e  = n/(n-2)*(w*(qe(2)-qe(1)))^2;
		% variance of the dependent variable
		qy  = quantile(Y,[p 1-p]);
		var_y = (w*(qy(2)-qy(1)))^2;
		% variance of the independent variable	
		qx  = quantile(Y,[p 1-p]);
		var_x = (w*(qx(2)-qx(1)))^2;
	end
	% standard error of the interecept
	s_inter = sqrt(var_e/n);
	% standard error of the slope
	s_slope = sqrt(var_e/(n*var_x));
	% coefficient of correlation
	R2 = 1 - var_e/var_y;
	
	c = [intercept; slope];
end % robustlinreg

