% Sat 18 Apr 22:02:13 +08 2020
%
% solve first order ordinary differntial equation as 
% quasi boundary value problem
% (intended to be used in combination with bvp2c for coubpled first and second
%  order equations)
function [x, y, out] = bvp1c(odefun,bcfun,xi,varargin)
	opt = bvp2_check_arguments(varargin{:});
	nx  = opt.nx;

	% segment end points
	x   = mesh1(xi,opt.nx,opt.xs);

	% segment mid points
	xc  = mid(x);

	% segment lengths
	dx  = diff(x);

	% number of segments
	nxc = nx-1;

	% number of equations per segment
	m   = 2;

	out = struct();

	% allocate memory for differential operator
	A   = sparse([],[],[],m*nxc,m*nxc,6*m*nxc);

	% allocate memory for inhomogeneous part
	b   = zeros(m*nxc,1);

	% get number of coupled odes
%	cc  = feval(odefun);
%	neq = size(cc,3);
	neq = 1;

	% initial value of solution at segment midpoints
	ypm    = rand(m*nxc*neq,1);
	ll     = [];

	% solve non-linear system by picard iteration
	% TODO gauss newton
	[ypm, out.cflag] = picard(@bvp1c_,ypm,opt.sopt);

	% interpolate solution to segment end points (grid points)
	y = inner2outer_bvp1c(ypm,ll);
	%y = inner2outer(ypm(1:2:end-1));
	%y = y + [ypm(2:2:end);ypm(end)];


function ypm = bvp1c_(ypm)

	% ode-coefficients at section mid-points
	cc = feval(odefun,xc(1),0);

	% function values at section mid-points
	yc = zeros(neq*nxc,1);	
	for idx=1:neq
		c = cc(1,end-2:end,idx);
		% TODO this silently assumes globally  degeneration / no degeneration
		if (0 ~= c(1,2,neq)) 
		yc((idx-1)*nxc+1:idx*nxc) = ( ...
			   		  ypm(2*nxc*(idx-1)+1:2:2*nxc*idx) ...
					+ ypm(2*nxc*(idx-1)+2:2:2*nxc*idx) );
		else
			yc((idx-1)*nxc+1:idx*nxc) = ... 
			   		ypm(2*nxc*(idx-1)+2:2:2*nxc*idx);
		end
	end

	% ode-coefficients at section mid-points
	cc = feval(odefun,xc,yc);


	% TODO allocate AA
	AA = [];
	rr = zeros(neq*m*nxc,1);
        %ii = zeros(neq*m*nxc,1);
	ll = zeros(nxc,1,neq);

	% for each of the coupled odes
	for ccdx=1:neq
		ll(:,:,ccdx) = roots1(cc(:,end-2:end-1,ccdx));

		[A, b] = bvp1c_assemble(cc,ll,ccdx,dx,xi,bcfun);

		% stack system of odes
		AA(1+(ccdx-1)*m*nxc:ccdx*m*nxc, ...
		   1+(ccdx-1)*m*nxc:ccdx*m*nxc) = A;
		rr((ccdx-1)*m*nxc+1:ccdx*m*nxc) = b;
		%ii((ccdx-1)*m*nxc+1:ccdx*m*nxc) = ih;
		%ii = 0;
	end % for idx
	% solve
	ypm = (AA \ rr); %- ii;
end

function y = inner2outer_bvp1c(ypm,ll)
	y = zeros((nxc+1)*neq,1);
	nx = nxc+1;
	for ccdx=1:neq
		r = ll(:,1,ccdx);
		if (0 ~= r(1,1,ccdx))
		y(m*(ccdx-1)*nx+1:ccdx*nx-1) = ( ...
		   myexp(-0.5*r(:,1,ccdx).*dx).*ypm(m*(ccdx-1)*nxc+1:2:m*ccdx*nxc) ...
		       + ypm(m*(ccdx-1)*nxc+2:2:m*ccdx*nxc) ...
			);
		y(ccdx*nx)                 = ( ...
			 myexp(+0.5*r(end,1,ccdx))*dx(end)*ypm(m*ccdx*nxc-1) ...
			 + ypm(m*ccdx*nxc) );
		else
			y_ = ypm(m*(ccdx-1)*nxc+2:2:m*ccdx*nxc);
			y((ccdx-1)*nx+1:ccdx*nx) = inner2outer(y_);
		end			      
	end
end

end % bvp2c_

