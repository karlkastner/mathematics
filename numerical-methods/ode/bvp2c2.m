% Mon 11 Mar 11:38:36 CET 2019
% Karl Kastner, Berlin
%
%% solve second order boundary value problem via roots of the characteristic
%% polynomial
%%
%% input:
%%
%% x   : [nx1] discretized domain
%%       n : number of vertices
%%      nxc = n-1 : number of segments
%%
%% bc  : struct : boundary condition
%%       bc.p(1)*y(0) + bc.pd(2)*y'(0) = bc.val(1)
%%       bc.p(2)*y(L) + bc.pd(2)*y'(L) = bc.val(2)
%%
%% output:
%%
%% A   : [2*nxc x 2*ns] disrcretisation matrix
%% rhs : [2*nxc x 1] right hand size
%%
%% y = A^-1 rhs
%%
% note: this can be simplifid by adding the two waves
%
function [x, y, out] = bvp2c2(odefun, bc, xi, opt)
	opt = bvp2_check_arguments(varargin{:});
	nx = opt.nx;

	% reach end points
	% TODO option for automatic mesh adaptation
	x = mesh1(xi,opt.nx,opt.xs);

	% reach mid points
	xc = mid(x);

	% reach lengths
	dx = diff(x);

	% number of segments
	nxc = nx-1;
	
	% number of dimenxcions	
	% TODO get from input
	m = 1;

	% allocate memory
	rhs  = zeros(2*m*nxc,1);
	A    = sparse([],[],[],2*m*nxc,2*m*ns,8*ns);
	E    = [];
	V    = [];
	Vi = [];

	% inial value for ypm
	ypm = zeros(2*m*nxc,1);
	y   = zeros(m*nxc,1);
	% solve non-linear system by picard iteration
	[ypm, out.cflag] = picard(@bvp2c2_,ypm,opt.sopt);

	ypm_ = reshape(ypm,2,[]).';
	% TODO compute derivative at grid points, if required

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
			%dx1 = x(k)   - xr(idx);
			%dx2 = x(k+1) - xr(idx);
			dx   = xr(idx) - xc(k);
			%yr(idx) =   ypm(3*(k-1)+1)*exp(dx*l(k,1)) ...
		        %          + ypm(3*(k-1)+3)*exp(dx*l(k,2));
			e = zeros(2,2);
			e(1,1) = exp(dx*squeeze(E(1,k)));
			e(2,2) = exp(dx*squeeze(E(2,k)));
			y_     = (V(:,:,k)*e*Vi(:,:,k))*ypm(2*k-1:2*k);
			yr(idx) = y_(1)+y_(2);
		end % for idx
	end % if opt.xr

function ypm = bvp2c2_(ypm)

	% ode coefficients at grid points
	cc  = feval(odefun,x,y);


	% for each dimenxcion
	for ccdx=1:neq
	odec = cc(:,:,ccdx);

	% roots of characteristic polynomials at each vertex
	%r = roots2(cc);
	% average
	% TODO evaluate at mid-point (or better integrate)
	rbar = 0.5*(r(1:end-1,:)+r(2:end,:));
	rbar = roots2(odec(:,1:3));

	% derivative
	dr_dx     = bsxfun(@times,diff(r,[],1),1./dx); 

	% relative derivative
	dr_dx_rel = bsxfun(@times,1./(rbar(:,2)-rbar(:,1)),dr_dx);
	% dr_dx_rel = 0*dr_dx_rel;

	% wave number
	% dy1/dx = T1 y1 + R2 y2
	% dy2/dx = R1 y1 + T2 y2
	%A = zeros(2,2,n-1);
	% TODO pass individually to eig
	A = [];
if (0)
	A(1,1,:) = -(-rbar(:,1) - dr_dx_rel(:,1));
	A(1,2,:) = -dr_dx_rel(:,2); 
	A(2,1,:) =  dr_dx_rel(:,1);
	A(2,2,:) =  (rbar(:,2) - dr_dx_rel(:,2));
end

	A(1,1,:) =  rbar(:,1) - dr_dx_rel(:,1);
	A(1,2,:) = -dr_dx_rel(:,2);
	A(2,1,:) =  dr_dx_rel(:,1);
	A(2,2,:) =  rbar(:,2) + dr_dx_rel(:,2);

	% eigenvalue decomposition
	[V, E] = eig2x2(A);
	Vi     = inv2x2(V);
	%E(1,:)    = -E(1,:);
	
	% continuity of value
	% value at -dx/2
	% yl = V*exp(-dx/2*E)*Vi;
	el        = zeros(2,2,nx-1);
	el(1,1,:) = exp(-0.5*dx.'.*squeeze(E(1,:)));
	el(2,2,:) = exp(-0.5*dx.'.*squeeze(E(2,:)));
	yl        = mtimes2x2(V,mtimes2x2(el,Vi));

	% continuity of derivative
	% derivative at -dx/2
	% dy_dxl = V*(E*exp(-dx/2)*E)*V;
	del        = zeros(2,2,nx-1);
	del(1,1,:) = E(1,:).*squeeze(el(1,1,:)).';
	del(2,2,:) = E(2,:).*squeeze(el(2,2,:)).';
	dyl_dx     = mtimes2x2(V,mtimes2x2(del,Vi));
	%scale     = ones(nxc,1);
	%1./mean(abs(E));
	%scale = 1./abs(E(1,:));
	scale = ones(nxc,1);

	% value at +dx/2
	er        = zeros(2,2,nx-1);
	er(1,1,:) = exp(+0.5*dx.'.*squeeze(E(1,:)));
	er(2,2,:) = exp(+0.5*dx.'.*squeeze(E(2,:)));
	yr        = mtimes2x2(V,mtimes2x2(er,Vi));
%er_ = er;
	% derivative at +dx/2
	der        = zeros(2,2,nx-1);
	der(1,1,:) = E(1,:).*squeeze(er(1,1,:)).';
	der(2,2,:) = E(2,:).*squeeze(er(2,2,:)).';
	dyr_dx     = mtimes2x2(V,mtimes2x2(der,Vi));

	rhs = zeros(2*nxc,1);
	A   = zeros(2*nxc,2*ns);
	
	% boundarcy condition at left end
	neq = 1;
	A(neq,1) =         bc(1).p(1,1,ccdx)*(yl(1,1,1) + yl(2,1,1)) ...
			 + bc(1).p(1,2,ccdx)*(dyl_dx(1,1,1) + dyl_dx(2,1,1));
	A(neq,2) =         bc(1).p(2,1,ccdx)*(yl(1,2,1) + yl(2,2,1)) ...
			 + bc(1).p(2,2,ccdx)*(dyl_dx(1,2,1) + dyl_dx(2,2,1));
	rhs(1)   =         bc(1).val(1,ccdx);

	% for each interior segment
	% 2*(nxc-1) = 2*ns-2 = 2*(np-2) = 2*np-4
	for id=1:nxc-1
		neq = neq+1;
		% y_i(+dxi/2) - y_i+1(-dx_i+1/2) == 0, y = yl+yr
		A(neq,2*id-1)     =  (yr(1,1,id)   + yr(2,1,id));
		A(neq,2*id)       =  (yr(1,2,id)   + yr(2,2,id));
		A(neq,2*(id+1)-1) = -(yl(1,1,id+1) + yl(2,1,id+1));
		A(neq,2*(id+1))   = -(yl(1,2,id+1) + yl(2,2,id+1));
		% rhs is zero

		% dy_i/dx(+dxi/2) - dy_{i+1}/dx(-dx_{i+1}/2) == 0
		neq = neq+1;
		A(neq,2*id-1)     =  scale(id)*(dyr_dx(1,1,id)   + dyr_dx(2,1,id));
		A(neq,2*id)       =  scale(id)*(dyr_dx(1,2,id)   + dyr_dx(2,2,id));
		A(neq,2*(id+1)-1) = -scale(id)*(dyl_dx(1,1,id+1) + dyl_dx(2,1,id+1));
		A(neq,2*(id+1))   = -scale(id)*(dyl_dx(1,2,id+1) + dyl_dx(2,2,id+1));
		% rhs is zero
	end % for id

	% boundary condition at right end
	neq           = neq+1;
	A(neq,2*nxc-1) =   bc(2).p(1,1,ccdx)*(yr(1,1,ns)     + yr(2,1,ns)) ...
			 + bc(2).p(1,2,ccdx)*(dyr_dx(1,1,nxc) + dyr_dx(2,1,ns));
	A(neq,2*nxc)   =   bc(2).p(2,1,ccdx)*(yr(1,2,ns)     + yr(2,2,ns)) ...
			 + bc(2).p(2,2,ccdx)*(dyr_dx(1,2,nxc) + dyr_dx(2,2,ns));
	rhs(neq)      =   bc(2).val(1,ccdx);

	% stack system of ODEs
	% TODO
	AA = A;
	rr = rhs;

	end % for ccdx

	% balance
	s  = 1; %1./abs(diag(AA));
	AA = diag(s)*AA;
	rr = s.*rr;

	% solve
	ypm = AA \ rr;


	% function value at grid points
	y = [   squeeze(el(1,1,:)).*ypm(1:2:end-1) ...
              + squeeze(el(2,2,:)).*ypm(2:2:end);
		squeeze(er(1,1,end)).*ypm(end-1) ...
	      + squeeze(er(2,2,end)).*ypm(end) ];
    end % bvp2c2_
end % bvp2c2

