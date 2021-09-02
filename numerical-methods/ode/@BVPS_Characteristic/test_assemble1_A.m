%else
	A   = sparse([],[],[],m*nxc,m*nxc,2*m*nxc);
	b   = zeros(m*nxc,1);
	for id=1:nxc-1
		% homogeneous part, continuity between right end of segment i
		% and left end of segment j
		if (0 ~= c(id,2))
			A(2*id-1,2*id-1)   = -obj.exp(+0.5*r(id).*dx(id));
			A(2*id-1,2*id)     = -1;
			A(2*id-1,2*id+1)   = +obj.exp(-0.5*r(id+1).*dx(id+1));
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
			A(2*id,2*id-1)     = 1;
			b(2*id)            = c(id,3)/c(id,1);
		end
		% rhs
	end
%end

