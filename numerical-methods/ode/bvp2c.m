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
% TODO option to approximate exp(x) as 1+x
function [x, y, out] = bvp2c(odefun,bcfun,xi,varargin)
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
		yc      = feval(opt.ifun,x);
		for idx = 1:neq
			% assign initial values into the inhomogeneous part
			% this does not matter, as yc is reassembled later
			% the inhomogeneous part is always the second part
			ypm(nxc*(mi(idx)+1)+1:nxc*(mi(idx)+2)) = yc((idx-1)*nxc+1:idx*nxc);
		end
	end
	ll     = [];

	% solve non-linear system by picard iteration
	[ypm, out.cflag, out.kiter] = picard(@bvp2c_,ypm,opt.sopt);

	% interpolate solution to segment end points (grid points)
	y = inner2out_(ypm,ll);

	% derivative of solution
	if (nargout > 3)
	dydx    = [ (  l(:,1).*ypm(1:m:m*nxc-m+1).*exp(-l(:,1).*dx(1:nxc)/2) ...
                     + l(:,2).*ypm(3:m:m*nxc    ).*exp(-l(:,2).*dx(1:nxc)/2));
	            (  l(nxc,1).*ypm(m*nxc-2)*exp(l(nxc,1)*dx(nxc)/2) ...
                     + l(nxc,2).*ypm(  m*nxc)*exp(l(nxc,2)*dx(nxc)/2)) ];
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
			yr(idx) =   ypm(3*(k-1)+1)*exp(dx*l(k,1)) ...
		                  + ypm(3*(k-1)+3)*exp(dx*l(k,2));
		end % for idx
	end % if opt.xr

function ypm = bvp2c_(ypm)
	% ode-coefficients at section mid-points
	cc = feval(odefun,xc(1),0);

	% function values at section mid-points
	yc = zeros(neq*nxc,1);
	for idx=1:neq
		% yc = ypm(1:m:end-m+1) + ypm(3:m:end);
		if (2 == mm(idx))
			% first order ode
			c = cc(1,end-2:end,idx);
			if (0 ~= c(1,2)) 
			% not degenerated
			yc((idx-1)*nxc+1:idx*nxc) = ( ...
				  ypm(mi(idx)*nxc+1:mm(idx):mi(idx+1)*nxc) ...
				+ ypm(mi(idx)*nxc+2:mm(idx):mi(idx+1)*nxc) );
			else
			yc((idx-1)*nxc+1:idx*nxc) = ... 
			   		ypm(mi(idx)*nxc+2:2:mi(idx+1)*nxc);
				
			end
		else
			% 2nd order ode
			% sum of left and right going wave as well as inhomogeneous part
			yc((idx-1)*nxc+1:idx*nxc) = ( ...
				  ypm(mi(idx)*nxc+1:mm(idx):mi(idx+1)*nxc) ...
				+ ypm(mi(idx)*nxc+2:mm(idx):mi(idx+1)*nxc) ...
				+ ypm(mi(idx)*nxc+3:mm(idx):mi(idx+1)*nxc) );
		end
	end

%	if (m > 2)
		% inhomogeneous part
		%yi = yc(m:m:end);
%	end
	
	% ode-coefficients at section mid-points
	cc = feval(odefun,xc,yc);

	% TODO allocate AA
	AA = [];
	rr = zeros(mi(end)*nxc,1);
	%ii = []; %zeros(neq*m*nxc,1);
	ll = zeros(nxc,2,neq);

	% for each of the coupled odes
	for ccdx=1:neq
		% eigenvalues, roots of the characteristic polynomial
		% roots are in general not complex conjugate pairs,
		% as this is only the case when coefficients are real
		% and case the solution a damped wave
		l            = roots2(cc(:,1:3,ccdx));
		ll(:,:,ccdx) = l;

		% silent asumption : ode are either purely 1st or 2nd order, but
		% not mixed
		if (0 == cc(1,1,ccdx))
			[A, b] = bvp1c_assemble(cc,ll,ccdx,dx,xi,bcfun);
		else
			[A, b] = bvp2c_assemble(cc,ll,ccdx,dx,xi,bcfun);
		end
	
		% stack system of odes
		% the equations are only weekly non-linear,
		% so they can be solved individually
		AA(1+mi(ccdx)*nxc:mi(ccdx+1)*nxc, ...
		   1+mi(ccdx)*nxc:mi(ccdx+1)*nxc) = A;
		rr(1+mi(ccdx)*nxc:mi(ccdx+1)*nxc,1) = b;
%		AA(1+(ccdx-1)*mc*nxc:ccdx*mc*nxc, ...
%		   1+(ccdx-1)*mc*nxc:ccdx*mc*nxc) = A;
%		rr((ccdx-1)*mc*nxc+1:ccdx*mc*nxc,1) = b;
%		ii((ccdx-1)*mc*nxc+1:ccdx*mc*nxc,1) = ih;
	end % for ccdx
%cA = condest(AA)
%pause

	% balance
	s  = 1; %1./abs(diag(AA));
	AA = diag(s)*AA;
	rr = s.*rr;
	% solve
	ypm    = (AA \ rr);

end % bvp2c_

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
		y(n0_+1:n0_+nx) = [(  ypm(n0+1:m:n0+m*nxc-m+1).*exp(-0.5*l(:,1).*dx(1:nxc)) ...
			              + ypm(n0+2:m:n0+m*nxc-m+2) ...
	                              + ypm(n0+3:m:n0+m*nxc).*exp(-0.5*l(:,2).*dx(1:nxc)) );
		                   (  ypm(n0+m*nxc-m+1)*exp(0.5*l(nxc,1)*dx(nxc)) ...
			              + ypm(n0+m*nxc-m+2) ...
	                              + ypm(n0+m*nxc)*exp(0.5*l(nxc,2)*dx(nxc))) ...
			];
		end
	end
	

function [A,b,ih] = bvp2c_assemble(cc,ll,ccdx,dx,xi,bcfun)
	odec = cc(:,:,ccdx);
	l    = ll(:,:,ccdx);
	m    = 3;
%	odec = bsxfun(@times,odec,1./odec(:,1));
	
	b = zeros(nxc,1);
	A = sparse([],[],[],m*nxc,m*nxc,(12+2)*nxc);

	% boundary condition at left end
	% f0'(0) = i o z1
	% alternatively, value and derivative could be specified
	% in separate rows, however, than this becomes an ivp and
	% no bc on the right side has to be specified
%	nout = nargout(bcfun);
	% TODO check sign of imaginary part to decide which root is the left going
	nout = 3;
	if (nout < 3)
		[v, p] = bcfun(xi(1),[],ccdx);
		q = [1, 1];
	else
		[v, p, q] = bcfun(xi(1),[],ccdx);
	end
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end
	A(1,1:m) = [ q(1)*(p(1) + p(2)*l(1,1))*exp(-0.5*l(1,1)*dx(1)), ...
		     p(1), ...
		     q(2)*(p(1) + p(2)*l(1,2))*exp(-0.5*l(1,2)*dx(1)), ...
		   ];
	% TODO scale all rows by maximum or diagonal value
	scale    = 1; %1./max(abs(A(1,1:m)));
	A(1,1:m) = scale*A(1,1:m);

	% right hand size
	% b    = -c(1:2:end-1,4);, nope, as this is not the ode, but y1 = y2
	b(1) = scale*v;

	% ode in interior sections
	% at each interior section end point
	for k=1:nxc-1
		% continuity of value
		% f_k(k dl) = f_k+1(-k dl)
		% rhs term must be in the centre, as first an last row are overwritten
		% 3,6,9
		A(m*(k-1)+3, m*(k-1) + (1:2*m)) = [	exp( l(k,1)*dx(k)/2), ...
						        1, ...
						        exp( l(k,2)*dx(k)/2), ...
			                               -exp(-l(k+1,1)*dx(k+1)/2), ...
					               -1, ...
					               -exp(-l(k+1,2)*dx(k+1)/2) ...
						  ];
		% continuity of first derivative
		% f'_k(k dl) = f'k+1(-k dl)
		scale = 1; %1./abs(l(k,1));
		% 4,7,10
		A(  m*(k-1)+4, m*(k-1) + (1:2*m)) = scale*[ ...
					   l(k,1)*exp(l(k,1)*dx(k)/2), ...
				         0, ...
					   l(k,2)*exp(l(k,2)*dx(k)/2), ...
					-l(k+1,1)*exp(-l(k+1,1)*dx(k+1)/2), ...
					 0, ...
					-l(k+1,2)*exp(-l(k+1,2)*dx(k+1)/2) ...
					 ];
	end % for k

	% ode at each section mid-point
	% TODO this is only true when c0 != 0
	% general solution : al*exp(+l*x) + ar*(exp(+r*x)) + a0 = 0
	for k=1:nxc
		% 2,5,8
		A(m*(k-1)+2, m*(k-1) + 1) = odec(k,1)*l(k,1).^2 + odec(k,2)*l(k,1) + odec(k,3);
		A(m*(k-1)+2, m*(k-1) + 2) = odec(k,3);
		A(m*(k-1)+2, m*(k-1) + 3) = odec(k,1)*l(k,2).^2 + odec(k,2)*l(k,2) + odec(k,3);
	end

	% rhs
	if (m>2)
		% 2,5,8,...
		b(m-1:m:m*nxc-1) = -odec(:,4);
	end

	% boundary condition at right end
	% fn(L)  = 0 (or better asymptotic y' = r y)
	if (nout < 3)
		[v, p] = bcfun(xi(2),[],ccdx);
		q      = [1, 1];
	else
		[v, p, q] = bcfun(xi(2),[],ccdx);
	end
	if ( 0 == sum(abs(p(1:2))) )
		error('weights must be non-zero');
	end
	A(m*nxc,m*nxc+(-m+1:0)) = [ ...
			q(1)*(p(1) + p(2)*l(nxc,1))*exp(+0.5*l(nxc,1)*dx(nxc)), ...
			p(1), ...
		        q(2)*(p(1) + p(2)*l(nxc,2))*exp(+0.5*l(nxc,2)*dx(nxc)), ...
			      ];
	
	% rhs
	b(m*nxc) = v;
	ih = 0;
end % bvp2c_assemble

end % bvp2c
