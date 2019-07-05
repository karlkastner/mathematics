% Sun 21 Aug 16:28:36 CEST 2016
% Karl Kastner, Berlin
%
%% numerical hessian projected to one dimenstion
%
function [Hp gp fc] = hessian_projected(fun,x0,dir,h,mode)
%	parflag = false; %true;
	if (nargin() < 4 || isempty(h))
		% note: h < eps^0.5 for consistent hessian !!!
		% h for gradient can be eps^0.5
		h = max(abs(x0)*eps.^0.25,eps.^0.25);
	end
	% make sure the direction is normalised (paranoid)
	dir = dir/norm(dir);

	% function value at origin
	fc = fun(x0);

	% left function value
	fl = fun(x0-dir*h);

	% right function value
	fr = fun(x0+dir*h);
	
	% projected gradient
	gp = 0.5*(fr - fl)/h;

	% porjected hessian
	Hp = (fl-2*fc+fr)/(h*h);
end

