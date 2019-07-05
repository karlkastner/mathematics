%Mon Nov 24 14:58:24 CET 2014
% Karl Kastner, Berlin
%% fit a power law
%% like polyfit, but returns parameter error estimates
% TODO parameter uncertainty
% TODO weighing in linearised form
function [param res A obj] = fit(obj,X,Y,W)

	if (nargin() < 4)
		W = 1;
	end
	X = double(X);
	Y = double(Y);
	W = double(W);

	X = cvec(X);

	% regression matrix
	A = [ones(size(X)) log(X)];

	% linear fit for initial parameters
	c_lin = A \ log(Y);
	
	% nonlinear fit
	sqrtW = sqrt(W);
	if (~obj.nonlinear)
	c = c_lin;
	else
	c = NaN(size(c_lin));
	for idx=1:size(c_lin,2)
		try
			c(:,idx) = lsqnonlin(@(c) sqrtW.*(exp(A*c) - Y(:,idx)),c_lin(:,idx));
		catch e
			e	
		end
	end
	end

	% residual
	res = exp(A*c) - Y;

	% root mean square error
	rmse = rms(res);

	% goodness of fit
	r2 = 1 - rmse.^2/var(Y);

	obj.param   = c;
	obj.rmse    = rmse;	
	obj.r2      = r2;

	obj.params = zeros(size(c));
	% Hessian and gradient
	% TODO weighing
	% g = 2*A'*exp(2*A*c) - 2*A'*(exp(A*c).*Y)
	for idx=1:size(c,2)
		H  = 4*(diag(exp(A*c(:,idx)))*A)'*diag(exp(A*c(:,idx)))*A;
		% in case of OLS H = 2*A'A
		%AA = 1/2*H;
		%C0 = inv(AA);
		% TODO why is this 1/4 H and not 1/2 HH (!)
		C0 = inv(1/4*H);
		obj.params(:,idx)  = rmse(idx)*sqrt(diag(C0));
	end

	%obj.C       = C;
	%obj.C0	    = C0;
	%obj.params  = sqrt(diag(C0)*serr2);
	%obj.nsample = ne;
	%obj.ssr     = ssr;
end % fit

