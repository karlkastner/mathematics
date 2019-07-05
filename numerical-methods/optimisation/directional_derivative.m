% Sun 21 Aug 16:28:36 CEST 2016
% Karl Kastner, Berlin
%
%% directional (projected) derivative
%% d : derivative, highest first
%% p : series expansion around x0
function [d p] = directional_derivative(fun,x0,dir,h,order,mode)
%	parflag = false; %true;
	if (nargin() < 4 || isempty(h))
		% note: h < eps^0.5 for consistent hessian !!!
		% h for gradient can be eps^0.5
		e = eps^(2^(1-order));
		h = max(norm(x0)*e,e);
		h_ = h;
	end

	% make sure the direction is normalised (paranoid)
	dir = dir/norm(dir);

	% number of support points
	n = double(idivide(int32(order+1),2));
	m = 2*n+1;

	% evaluate function at lagrangian support points
	f = zeros(m,1);
	h = h*(-n:n)';
	for idx=1:m
		f(idx) = feval(fun,x0 + dir*h(idx));
	end

	% vandermonde matrix
	A = vander_1d(h,order);

	% polynomial coefficients
	p = A \ f;

	p = flipud(p);

	% derive and evaluate polynomials n times and evaluate at zero

	d = p.*factorial((order:-1:0))';

if (0)
	d
	h=h_
	d_=[[0 0 1 0 0]*f;
         1*[1 -8 0 8 -1]./(12*h)*f
	 2*[-1 16 -30 16 -1]./(24*h^2)*f
	 6*[-1 2 0 -2 1]./(12*h^3)*f
	 24*[1 -4 6 -4 1]./(24*h^4)*f];
	flipud(d_)
pause
end
%	% function value at origin
%	fc = fun(x0);
%
%	% left function value
%	fl = fun(x0-dir*h);
%
%	% right function value
%	fr = fun(x0+dir*h);
%	
%	% projected gradient
%	gp = 0.5*(fr - fl)/h;
%
%	% porjected hessian
%	Hp = (fl-2*fc+fr)/(h*h);
end

