% Sun  8 Oct 11:05:23 CEST 2017
% Karl Kastner, Berlin
% 
%% solve system of non-linear second order odes (in more than one variable)
%% as boundary value problems by the finite difference method
%%
%% odefun provides ode coefficients c:
%% c(x,1) y''(x) + c(x,2) y'(x) + c(x,3) y = c(x,4)
%%  c_1 y" + c_2 y' + c_3 y + c_4 = 0
%%
%% subject to the boundary conditions
%% bcfun provides v and p and optionally q, so that:
%%
%% b_1 y + b_2 y' = f
%%    q(x,1)*( p(x,1) y_l(x) + p(x,2)  y_l'(x)
%%  + q(x,2)*( p(x,1) y_r(x) + p(x,2) y_r'(x)    = v(x)
%% where q weighs the waves travelling from left to right and right to left (default [1 1])
%
% TODO the system of odes is not yet coupled
function [x, y, cflag] = bvp2fdm(odefun,bcfun,X,opt)
	if (nargin()<4)
		opt = struct();
	end
	if (~isfield(opt,'nx'))
		nx = 100;
	else
		nx = opt.nx;
	end
	if (~isfield(opt,'xs'))
		xs = 1;
	else
		xs = opt.xs;
%		if (xs ~= 1)
%			warning('only equispaced grids are supported')
%			xs = 1;
%		end
	end
	if (~isfield(opt,'sopt'))
		sopt=struct();
	else
		sopt = opt.sopt;
	end
	if (~isfield(sopt,'maxiter'))
		sopt.maxiter = nx;
	end
	if (~isfield(sopt,'relaxation'))
		sopt.relaxation = 0.5;
	end
	if (~isfield(opt,'opc'))
		% second order accurate implementation of boundary conditions
		opt.obc = 1;
	end
	nx = max(nx,2);

	if (0)
		L  = X(2)-X(1);
		dx = L/(nx-1);
		x  = X(1) + dx*(0:nx-1)';
	end
	x  = mesh1(X,nx,xs);
	% step length
	dx = diff(x);

	% difference matrices
if (1)
	D2  = derivative_matrix_2_1d(x,[]);
	D1c = derivative_matrix_1_1d(x,[],2);
	D1r = derivative_matrix_1_1d(x,[],'+1');
	D1l = derivative_matrix_1_1d(x,[],'-1');
else
	D2  = derivative_matrix_2_1d(nx,L);
	D1c = derivative_matrix_1_1d(nx,L,2);
	D1l = derivative_matrix_1_1d(nx,L,'+2');
	D1r = derivative_matrix_1_1d(nx,L,'-2');
end

if (0)
	D2  = 1/dx^2*spdiags(ones(nx,1)*[ 1 -2 1],-1:1,nx,nx);
	D1c = 0.5/dx*spdiags(ones(nx,1)*[-1  0 1],-1:1,nx,nx);
	D1l = 1/dx*spdiags(ones(nx,1)*[-1  1],-1:0,nx,nx);
	D1r = 1/dx*spdiags(ones(nx,1)*[-1  1],0:1,nx,nx);
	%full(D2-D2_)
%	full(D1l)
%	full(D1l_)
%	full(D1r-D1r_)
	%full(D1c-D1c_)
%	pause
end

	I   = speye(nx);

	% get number of coupled odes
	cc   = feval(odefun);
	neq = size(cc,3);

	% initial value of y
	if (~isfield(opt,'icfun'))
		y = zeros(nx*neq,1);
	else
		y = opt.icfun(x);
	end

	% solve non-linear system by picard iteration
	[y, cflag] = picard(@bvp2fdm_,y,sopt);

function y = bvp2fdm_(y)
	cc = feval(odefun,x,y);
	AA = [];
	bb = [];

	% for each dimension
	for ccdx=1:neq
	c = cc(:,:,ccdx);

	% ode discretisation matrix in interior of computational domain
	A =	  diag(sparse(c(:,1)))*D2 ...
		+ diag(sparse(c(:,3)));

	% for the first derivative (advection term)
	% switch to forward/backward differencing when c1 == 0 and c3 == 0
	cdx =   abs(c(:,1)) > sqrt(eps) ...
	      | abs(c(:,3)) > sqrt(eps);
	%cdx = ones(size(cdx));
	%r   = roots(c(:,2:3)).*c(:,4);
	ldx = ~(c(:,2).*c(:,4) > 0);
	%ldx = (~cdx).*(r<0);
	%ldx = (~cdx).*(r>0);
	A   = A + diag(sparse(c(:,2)))*( ...
		  diag(sparse(cdx))*D1c ...
		+ diag( sparse((~cdx).*(ldx)) )*D1l ...
		+ diag( sparse((~cdx).*(~ldx)) )*D1r );
%	(*D1 ...

	% rhs
	b = -c(:,4);	

	% apply robin boundary conditions
	% f(0)*p1 + p2*f'(0) = f0;
	[f, p]     = bcfun(X(1),y(1),ccdx);

	p(end+1:3) = 0;
	A(1,:)    = 0;
%p
%	% improve conditioning for D"
%	if (0 == p(3))
%		if (0 == p(2))
%			scale = 1/dx(1)^2;
%		else
%			scale = 1/dx(1);
%		end
%	else
%		scale = 1;
%	end
%	scale = 1./scale;

	% TODO: second order correction not working for variable grid spacing
	if (opt.obc>1 && nx > 2)
		% third order
		A(1,1:3)  = [p(1) - 3/2*p(2)/dx(1) +   p(3)/dx(1)^2, ...
                                      2*p(2)/dx(1) - 2*p(3)/dx(1)^2, ...
                                   -1/2*p(2)/dx(1) +   p(3)/dx(1)^2];
	else
		A(1,1:3)  = [p(1) - p(2)/dx(1) +   p(3)/dx(1)^2, ... 
                                    p(2)/dx(1) - 2*p(3)/dx(1)^2, ...
					    +   p(3)/dx(1)^2 ...
                            ];
	end
	b(1)      = f;
%	full(A(1:2,1:3))
%	full(A(1,1:3))*y(1:3)
%	b(1:2)
%'honk'
%pause

	% f(L)*p00 + p01*f'(L) = fL;
	[f, p]     = bcfun(X(2),y(end),ccdx);
%	if (0 == p(3))
%		if (0 == p(2))
%			scale = 1/dx(end)^2;
%		else
%			scale = 1/dx(end);
%		end
%	else
%		scale = 1;
%	end
%	scale = 1./scale;
	if (~isempty(f))
	p(end+1:3) = 0;
	A(end,:)  = 0;
	if (opt.obc > 1 && nx > 2)
		A(end,end-2:end) = [       1/2*p(2)/dx(end) +   p(3)/dx(end)^2, ...
				            -2*p(2)/dx(end) - 2*p(3)/dx(end)^2, ...
                                    p(1) + 3/2*p(2)/dx(end) +   p(3)/dx(end)^2];
	else
		A(end,end-2:end) = [                      p(3)/dx(end)^2, ...
				             -p(2)/dx(end) - 2*p(3)/dx(end)^2, ...
                                       p(1) + p(2)/dx(end) +   p(3)/dx(end)^2];
	end
	b(end)    = f; 
	end

	% stacksystem of odes
	AA(1+(ccdx-1)*nx:ccdx*nx,1+(ccdx-1)*nx:ccdx*nx) = A;
	bb   = [bb; b];

	end % for ccdx

	% solve
	y = AA\bb;

end % bvp2fdm_

end % bvp2fdm

