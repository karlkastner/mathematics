% Sat 18 Apr 22:02:13 +08 2020
function [A,b] = bvp1c_assemble(cc,ll,ccdx,dx,xi,bcfun)
	% only use last 3 coefficients
	c  = cc(:,end-2:end,ccdx); 
	m  = 2;

	% root for homogeneous part
	%r = -c(:,2)/c(:,1);
	r  = ll(:,1,ccdx);

	% inhomogeneous part signed for rhs
	%ih = c(:,3)./c(:,2);
	ih = [];

	% match at rightend at segment i with right end at segment i+1
	nxc = size(c,1);
	A   = sparse([],[],[],m*nxc,m*nxc,2*m*nxc);
	b   = zeros(m*nxc,1);
	for id=1:nxc-1
		% homogeneous part, continuity between right end of segment i
		% and left end of segment j
		if (0 ~= c(id,2))
			A(2*id-1,2*id-1)   = -exp(+0.5*r(id).*dx(id));
			A(2*id-1,2*id)     = -1;
			A(2*id-1,2*id+1)   = +exp(-0.5*r(id+1).*dx(id+1));
			A(2*id-1,2*id+2)   = +1;
			% inhomogeneous part
			A(2*id,2*id)       = 1;
			b(2*id)            = -c(id,3)./c(id,2);
		else
			% degenerated, without damping,
			% the piecewise solution is a linear function
			% match left with right
			A(2*id-1,2*id-1)   = -0.5*dx(id);
			A(2*id-1,2*id)     = -1;
			A(2*id-1,2*id+1)   = -0.5*dx(id+1);
			A(2*id-1,2*id+2)   = +1;
			% slope
			A(2*id,2*id-1)       = 1;
			b(2*id)              = c(id,3)/c(id,1);
		end
		% rhs
		%b(id)      = (ih(id+1)-ih(id));
	end

	% boundary condition at left
	% only one boundary has to be specified for first order odes
	[v, p] = bcfun(xi(1),[],ccdx);
	nb = 0;
	if (~isempty(v))
		% p*f(0) + (1-p)*f(0)'' = v
		if (0~=c(1,2))
			A(2*nxc-1,1) = p*exp(-0.5*r(1)*dx(1)) + (1-p)*r(1)*exp(-0.5*r(1)*dx(1));
			A(2*nxc-1,2) = p;
			b(2*nxc-1)   = v;
		else
			A(2*nxc-1,1) = p*(-0.5*dx(1)) + (1-p);
			A(2*nxc-1,2) = p;
			b(2*nxc-1)   = v;
		end
		%b(2*nxc-1)   = p*ih(1) + v;
		nb           = 1;
	end

	[v, p] = bcfun(xi(2),[],ccdx);
	if (~isempty(v))
		if (0~=c(ndx,2))
			A(2*nxc-1,nxc-1) = (       p*exp(+0.5*r(nxc)*dx(nxc)) ...
					     + (1-p)*r(nxc)*exp(+0.5*r(nxc)*dx(nxc)) ...
					   );
			A(2*nxc-1,nxc)     = p;
			b(2*nxc-1)         = v;
		else
			A(2*nxc-1,nxc-1) = ( p*0.5*dx(nxc) + (1-p) );
			A(2*nxc-1,nxc)   = p;
			b(2*nxc-1)       = v;
		end
		nb     = nb+1;
	end
	if (0~=c(end,2))
		A(2*nxc,2*nxc) = 1;
		b(2*nxc)       = -c(nxc,3)./c(nxc,2);
	else
		% slope of last segment
		A(2*nxc,2*nxc-1) = 1;
		b(2*nxc)       = c(nxc,3)./c(nxc,1);
	end
	if (1 ~= nb)
		error('need boundary condition at exactly one end');
	end
%	full(A)
%	eig(full(A))
%AA = A(2:2:end,2:2:end);
%	full(AA)
%	eig(full(AA))
%pause
end % bvp1c_assemble

