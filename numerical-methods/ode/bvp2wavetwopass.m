% Mon  4 Dec 14:10:05 CET 2017
% Karl Kastner, Berlin
%
%% two pass solution for the linearised wave equation
%% solve first for the wave number k, and then for y
%
%function [x, z, q, k] = wave_twopassq(K2fun,Xi,omega,opt)
function [x, y, cflag, k] = bvp2wavetwopass(odefun,bcfun,Xi,opt)

	if (nargin()<4)
		opt = struct();
	end

	if (~isfield(opt,'nx'))
		nx = 1024;
	else
		nx = opt.nx;
	end

	% domain length
	L  = Xi(2)-Xi(1);
	dx = L/(nx-1);
	x  = Xi(1)+dx*(0:nx-1)';
	sopt.MaxStep = dx;
	
	%xy  = Xi(1)+dx*(0:nx-1)';
	y = zeros(nx,1);
	[y cflag] = picard(@wave_twopassq_,y);

	% resanple y and k

	% resample k to x
	%y = interp1(xy,y,x);
	if (nargout()>3)
		%k     = kfun(x);
		k     = interp1(xk,k,x);
	end

	% TODO, this is again not second order accurate
	% dq_dx = cdiff(q)./cdiff(x);
	% z     = -1./(1i*omega)*dq_dx;

function y = wave_twopassq_(y)
	flag = 0;

	% ode coefficients at L
	cL = odefun(Xi(2),y(end));
	% normalize coefficients
	cL = cL/cL(1);

	% roots of the characteristic polynomial
	% aka wave number in case of constant coefficients
	rL = roots2(cL);

	% initial value for k at right end (L)
	k_L = rL(2);
%	k_L^2 + cL(3)
%	k_L = sqrt(-cL(3))
%pause

	% note: non stiff solver leads to spurious oscillation
	solver_ = @ode23s;

	% from L to zero: solve for k
	[xk, k] = solver_(@kdot,[Xi(2),Xi(1)],k_L,sopt);
%	[xk_, k_] = ode23(@kdot,[Xi(2),Xi(1)],rL(2),opt);
%figure()
%subplot(2,2,1)
%kk= [k,interp1(xk_,k_,xk)];
%plot(xk,real(kk));
%subplot(2,2,2)
%plot(xk,imag(kk));
%pause

	% flip x back
	if (flag)
		xk = flipud(L-xk);
		k  = flipud(k);
	end

	% resample to equidistant grid
	% TODO use piecewise polynomial
	% k = interp1(xk,k,x);
	%dkdx  = cdiff(k)./cdiff(x);
	% kfun  = @(xi) interp1(x,k,xi);
	%dkfun = @(xi) interp1(x,dkdx,xi);

	% from 0 to L: solve for Q1
	Q10     = 1;
	[xy, y] = solver_(@ydot,Xi,Q10,sopt);

	% resample y
	%yfun    = @(x) interp1(xy,y,x);
	y = interp1(xy,y,x);

	% apply boundary condition
	[bv bp] = bcfun(1,y(1));
	scale   = bv/(bp(1)*y(1) + bp(2)*(y(2)-y(1))/dx);
	y       = scale*y;

	% apply bc
	%scale = 1/z(1);
	%z = scale*z;
	%q = scale*q;
end % wave_twopass
	
	function kdot = kdot(x_,k_)
		% solve backward in space
		if (flag)
			x_   = L-x_;
		end
		% ode-coefficients
		% TODO, preduild interpolators
		y_ = interp1(x,y,x_);
		c = odefun(x_,y_);
		% normalize
		c = c/c(1);
%		R  = roots2(c);
%		R2 = R(1).^2;
%		R2 = c(3);
%		c(3)
%		k_.^2
%pause
%		kdot = -(k_.^2 - R2);
		kdot = -(c(1)*k_.^2 + c(2)*k_ + c(3));
		%kdot = -(k_.^2 + R2);
	end % kdot
	
	% y is tidal discharge in case of river-tide
	function ydot = ydot(x_, y_)
		% qdot  = 1./kfun(x)*(K2fun(x)-dkfun(x)).*q;
		% ydot = kfun(x_).*y_;
		k_ = interp1(xk,k,x_);
		ydot = k_.*y_;
	end	
end % wave_twostep

