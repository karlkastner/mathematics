% Sat 18 Apr 22:02:13 +08 2020
function [A,b] = assemble1_A(obj, cdx, edx)
%cc,ll,ccdx,dx,xi,bcfun,bcarg)

	% only use last 3 coefficients
	c  = obj.out(cdx).cc(:,end-2:end,edx); 
	m  = 2;

	% root for homogeneous part
	r  = obj.out(cdx).ll(:,1,edx);

	% match at rightend at segment i with right end at segment i+1
	nxc = size(c,1);
if (1)
	nbuf = 0;
	Abuf = zeros(5*(nxc-1)+1,3);
	b    = zeros(m*nxc,1);

	id = (1:nxc-1)';
	% TODO if c==0 then r==0, exploit
	fdx = (c(:,2) ~= 0);
	Abuf(id,1) = 2*id-1;
	Abuf(id,2) = 2*id-1;
	% continuity between right end of segment i
	% and left end of segment j
	% degenerated, without damping,
	% the piecewise solution is a linear function
	% match left with right
	Abuf(id,3)   = (fdx(1:end-1).*-obj.exp(+0.5*r(id).*dx(id)) ...
			+ (1-fdx(1:end-1)).*-0.5.*dx(id));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+id,1) = 2*id-1;
	Abuf(nbuf+id,2) = 2*id;
	Abuf(nbuf+id,3) = -1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+id,1) = 2*id-1;
	Abuf(nbuf+id,2) = 2*id+1;
	Abuf(nbuf+id,3)  = (     fdx(2:end).*+obj.exp(-0.5*r(id+1).*dx(id+1)) ...
	                    + (1-fdx(2:end)).*-0.5.*dx(id+1));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+id,1) = 2*id-1;
	Abuf(nbuf+id,2) = 2*id+2;
	Abuf(nbuf+id,3) = +1;
	nbuf = nbuf+nxc-1;
	% inhomogeneous part
	id = (1:nxc)';
	Abuf(nbuf+id,1) = 2*id;
	Abuf(nbuf+id,2) = 2*id + (fdx-1);
	Abuf(nbuf+id,3) = 1;
	nbuf            = nbuf+nxc-1;
	b(2*id(fdx))    = -c(id(fdx),3)./c(id(fdx),2);
	b(2*id(~fdx))   = +c(id(~fdx),3)./c(id(~fdx),1);

	A = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),m*nxc,m*nxc);

%	if (0~=c(end,2))
%		A(2*nxc,2*nxc) = 1;
%		b(2*nxc)       = -c(nxc,3)./c(nxc,2);
%	else
%		% slope of last segment
%		A(2*nxc,2*nxc-1) = 1;
%		b(2*nxc)       = c(nxc,3)./c(nxc,1);
%	end


else
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
end

	% boundary condition at left
	% only one boundary has to be specified for first order odes
	%[v, p, ~, set] = bcfun(xi(1),[],ccdx);
	[v, p, ~, set] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});

	nb = 0;
	if (set) %~isempty(v))
		% p*f(0) + (1-p)*f(0)'' = v
		row = 2*nxc-1;
		if (0 ~= c(1,2))
			A(ow,1)  = p*obj.exp(-0.5*r(1)*dx(1)) + (1-p)*r(1)*obj.exp(-0.5*r(1)*dx(1));
			A(row,2) = p;
			b(row)   = v;
		else
			A(row,1) = p*(-0.5*dx(1)) + (1-p);
			A(row,2) = p;
			b(row)   = v;
		end
		nb           = 1;
	end

	%[v, p, ~, set] = bcfun(xi(2),[],ccdx);
	[v, p, ~, set] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});

	if (set) %~isempty(v))
		row=2*nxc-1; % -1
		if (0~=c(end,2))
			A(row,2*nxc-1) = (       p*obj.exp(+0.5*r(nxc)*dx(nxc)) ...
					     + (1-p)*r(nxc)*obj.exp(+0.5*r(nxc)*dx(nxc)) ...
					   );
			A(row,2*nxc)     = p;
			b(row)         = v;
		else
			A(row,2*nxc-1) = ( p*0.5*dx(nxc) + (1-p) );
			A(row,2*nxc)   = p;
			b(row)       = v;
		end
		nb     = nb+1;
	end

	if (1 ~= nb)
		error('need boundary condition at exactly one end');
	end
end % assemble1_A
