% Sat 28 Oct 14:43:06 CEST 2017
% Karl Kastner, Berlin
%
%% solve system of non-linear second order odes (in more than one variable)
%% as boundary value problems
%%
%% odefun provides ode coefficients c:
%% c(x,1) y''(x) + c(x,2) y'(x) + c(x,3) y = c(x,4)
%%    c_1 y"     + c_2 y'       + c_3 y + c_4 = c_4
%%
%% subject to the boundary conditions
%% bcfun provides v and p and optionally q, so that:
%%
%% b_1 y + b_2 y' = f
%%    q(x,1)*( p(x,1) y_l(x) + p(x,2)  y_l'(x)
%%  + q(x,2)*( p(x,1) y_r(x) + p(x,2) y_r'(x)    = v(x)
%% where q weighs the waves travelling from left to right and right to left (default [1 1])
%
% TODO coupling for multi-reach odes (several x)
% TODO use buffers instead of sparse
% TODO option for automatic mesh adaptation
% TODO option to approximate myexp(x) as 1+x
function [x, y, out] = bvp2c(odefun,bcfun,xi,varargin)
	myexp = @(x) (1+x);

	opt = bvp2_check_arguments(varargin{:});
	nx = opt.nx;

	% segment end points
	x = mesh1(xi,opt.nx,opt.xs);

	% segment mid points
	xc = mid(x);

	% segment lengths
	dx  = diff(x);

	% number of segments
	nxc = nx-1;

	% number of equations per segment

	out = struct();
	
	% determine number of coupled odes
	cc  = feval(odefun,xi(1),1);
	neq = size(cc,3);
	% determine if ode is 1st (m=2) or 2nd order (m=3)
	mm  = 2*ones(neq,1) + (0 ~= flat(cc(1,1,:)));
	mi = cumsum([0;mm]);

	% allocate memory for differential operator
	A  = sparse([],[],[],mi(end)*nxc,mi(end)*nxc,4*mi(end)*nxc);
	% allocate memory for inhomogeneous part
	b  = zeros(mi(end)*nxc,1);

	% initial value of ypm
	% complex amplitude of the left and right going wave at segment mid points
	% ypm    = zeros(m*nxc*neq,1);
	% initial value
	% TODO should reassembly not take place inside iteration?
	ypm    = zeros(mi(end)*nxc,1);
	if (isfield(opt,'ifun'))
		yi      = feval(opt.ifun,x);
		for idx = 1:neq
			% assign initial values into the inhomogeneous part
			% this does not matter, as yi is reassembled later
			% the inhomogeneous part is always the second part
			yci = mid(yi((idx-1)*(nxc+1)+1:idx*(nxc+1)));
			ypm(nxc*mi(idx)+2:mm(idx):nxc*mi(idx+1)) = yci;
		end
	end
	ll     = [];

	% solve non-linear system by picard iteration
	[ypm, out.cflag, out.kiter] = picard(@bvp2c_solve,ypm,opt.sopt);

	% interpolate solution to segment end points (grid points)
	y = inner2out_(ypm,ll);


	% derivative of solution
	if (nargout > 3)
	% TODO, this is not working any more, due to regression
	% since first and second order odes were coupled
	l = [];
	m = 0;
	dydx    = [ (  l(:,1).*ypm(1:m:m*nxc-m+1).*myexp(-l(:,1).*dx(1:nxc)/2) ...
                     + l(:,2).*ypm(3:m:m*nxc    ).*myexp(-l(:,2).*dx(1:nxc)/2));
	            (  l(nxc,1).*ypm(m*nxc-2)*myexp(l(nxc,1)*dx(nxc)/2) ...
                     + l(nxc,2).*ypm(  m*nxc)*myexp(l(nxc,2)*dx(nxc)/2)) ];
	end

	% resampling of solution
	if (isfield(opt,'xr') && ~isempty(opt.xr))
		xr = opt.xr;
		k  = 1;
		yr = zeros(size(xr));
		for idx=1:length(xr)
			while(xr(idx) > x(k+1) && k+1<nx)
				k = k+1;
			end % while
			% expand
			dx   = xr(idx) - xc(k);
			yr(idx) =   ypm(3*(k-1)+1)*myexp(dx*l(k,1)) ...
		                  + ypm(3*(k-1)+3)*myexp(dx*l(k,2));
		end % for idx
	end % if opt.xr



function ypm = bvp2c_solve(ypm)
	[AA,rr,ll] = bvp2c_assemble(ypm,odefun,bcfun,xi,xc,dx,neq,nxc,mm,mi);

	% balance
	s  = 1; %1./abs(diag(AA));
	AA = diag(s)*AA;
	rr = s.*rr;
	% solve
	ypm    = (AA \ rr);
end % bvp2c_solve


	function y = inner2out_(ypm,ll)
		y = zeros(nx*neq,1);
		for id=1:neq
			r = ll(:,1,id);
		     %n0  = m*nxc*(ne-1);
		     %n0_ = nx*(ne-1);
		     %y(n0_+1:n0_+nx) = inner2outer(ypm(n0+1:m:n0+m*nxc-m+1));
		     if (2 == mm(id))
			if (0 ~= r(1,1,id))
			y_ = (   ypm(mi(id)*nxc+1:mm(id):mi(id+1)*nxc) ...
			       + ypm(mi(id)*nxc+2:mm(id):mi(id+1)*nxc) );
			else
			% degenerated linear function
			y_ = ypm(mi(id)*nxc+2:mm(id):mi(id+1)*nxc);
			end
			y(nx*(id-1)+1:nx*id) = inner2outer(y_);
		     else
			y_ = (   ypm(mi(id)*nxc+1:mm(id):mi(id+1)*nxc) ...
			       + ypm(mi(id)*nxc+2:mm(id):mi(id+1)*nxc) ...
			       + ypm(mi(id)*nxc+3:mm(id):mi(id+1)*nxc) ...
		             );
			y(nx*(id-1)+1:nx*id) = inner2outer(y_);

		     %y(nx*(id-1)+1:nx*id) = ( ...
			%	  inner2outer(ypm(mi(id)*nxc+1:mm(id):mi(id+1))) ...
			%	+ inner2outer(ypm(mi(id)*nxc+2:mm(id):mi(id+1))) ...
			%	+             ypm(mi(id)*nxc+3:mm(id):mi(id+1)) ...
			%	);
		     end
		end
	end
	function y = inner2outer_bvp2c(ypm)
		y = zeros(nx*neq,1);
		for ne=1:neq
		     n0  = m*nxc*(ne-1);
		     n0_ = nx*(ne-1);
		     l = ll(:,:,ne);
		     y(n0_+1:n0_+nx) = [(  ypm(n0+1:m:n0+m*nxc-m+1).*myexp(-0.5*l(:,1).*dx(1:nxc)) ...
		  	                   + ypm(n0+2:m:n0+m*nxc-m+2) ...
	                                   + ypm(n0+3:m:n0+m*nxc).*myexp(-0.5*l(:,2).*dx(1:nxc)) );
		                        (  ypm(n0+m*nxc-m+1)*myexp(0.5*l(nxc,1)*dx(nxc)) ...
			                   + ypm(n0+m*nxc-m+2) ...
	                                   + ypm(n0+m*nxc)*myexp(0.5*l(nxc,2)*dx(nxc))) ...
			               ];
		end
	end

end % bvp2c

