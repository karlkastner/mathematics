%Mo 3. Aug 10:15:17 CEST 2015
% Karl Kastner, Berlin
%% draw random variables of two corrlated ar1 processes
function [x y]=randar1_dual(s1,s2,rho1,rho2,rho12,n,m)
	if (nargin() < 7)
		m = 1;
	end

	ex = randn(n,m);
	ey = randn(n,m);

	% scale factor of rho12
	s     = sqrt(1-rho1^2)*sqrt(1-rho2^2)/(1-rho1*rho2);
	rho12 = rho12/s;
	if (rho12 > 1)
		% Note that for high correlation rho12 larger one and the series becomes complex
		warning('condition cannot be satisfied');
	end
	A = [1,           rho12;
             0, sqrt(1-rho12^2)];
	
	% cross-correlate residuals
	% e_corr = A*[ex ey];
	ex_ = ex*A(1,1) + ey*A(2,1);
	ey_ = ex*A(1,2) + ey*A(2,2);

	% initial condition x(0) = 0, y(0) = 0
	ex_(1) = ex_(1)/sqrt(1-rho1^2);
	ey_(1) = ey_(1)/sqrt(1-rho2^2);

	% auto-correlate residuals
	% x(idx,:) = rho1*x(idx-1,:) + sqrt(1-rho1^2)*ex(idx,:);
	% y(idx,:) = rho2*y(idx-1,:) + sqrt(1-rho2^2)*ey(idx,:);
	x = filter(s1*sqrt(1-rho1^2),[1 -rho1],ex_);
	y = filter(s2*sqrt(1-rho2^2),[1 -rho2],ey_);
end % randar1_dual

