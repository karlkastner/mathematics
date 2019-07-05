% Fri 15 Mar 20:03:28 CET 2019
%
%% transformation matrix of second order ode
%% to left and right going wave
%%
%% c = odefun(x)
%% c1 y'' + c2' y + c3 y == 0
%% y = y_p + y_m, left and right going wave
%% d/dx [y_p, y_m] = A*[y_m, y_p]
%
% function [A, Ac, r, dr_dx] = ode2_matrix(odefun,x,y,dx,irule)
function [A, Ac, r, dr_dx] = ode2_matrix(odefun,x,y,dx,irule)
	if (2 == length(y))
		y = y(1)+y(2);
	end
	if (nargin()<3||isempty(dx)) % TODO use grad
		dx  = 1e-3; %*(x(end)-x(;
	end
	if (isempty(irule))
		irule = 't';
	end
	c_r = odefun(x+dx/2,y);
	r_r = roots2(c_r);
	%r_r = -(r_r);
%	r_r = [r_r(2),r_r(1)];
	c_l = odefun(x-dx/2,y);
	r_l = roots2(c_l);
	%r_l = -(r_l);
%	r_l = [r_l(2),r_l(1)];
	switch (lower(irule(1)))
%		case {'l'} % left
%			r = r_l;
%		case {'r'} % right
%			r = r_r;
%		case {'m'} % midpoint
%			c = odefun(x,y);
%			r = roots2(c);	
		case {'t'} % trapezoidal
			r = 0.5*(r_l + r_r);
		otherwise
			error('here');
	end

	dr_dx     = (r_r-r_l)/dx;
	%dr_dx_rel = 1/(r(2)-r(1))*dr_dx;
	dr_dx_rel = 1/(r(1)-r(2))*dr_dx;
	%dr_dx_rel = -conj(dr_dx_rel);
	%dr_dx_rel = -dr_dx_rel;

%	A = [r(1) + dr_dx_rel(1),        +dr_dx_rel(2);
%                 - -dr_dx_rel(2), r(2) - -dr_dx_rel(1)];
	A = [r(1) - dr_dx_rel(1),      - dr_dx_rel(2);
                  + dr_dx_rel(1), r(2) + dr_dx_rel(2)];
	Ac = [r(1),0;
              0,r(2)];
%	A = [r(1), 0;
%	     0   , -r(2)];
end

