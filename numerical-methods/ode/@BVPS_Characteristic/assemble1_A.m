% Sat 18 Apr 22:02:13 +08 2020
function assemble1_A(obj, cdx, edx)

	% only use last 3 coefficients
	c  = obj.out(cdx).cc(:,end-2:end,edx); 
	m  = 2;

	% root for homogeneous part
	r  = obj.out(cdx).ll(:,1,edx);

	% match at rightend at segment i with right end at segment i+1
	nxc = size(c,1);

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
	id = (1:nxc)'; % 2 : 2*nxc
	Abuf(nbuf+id,1) = 2*id;
	Abuf(nbuf+id,2) = 2*id + (fdx-1);
	Abuf(nbuf+id,3) = 1;
	% last to zero
	Abuf(nbuf+id(end),3) = 0;
	nbuf            = nbuf+nxc;
	b(2*id(fdx))    = -c(id(fdx),3)./c(id(fdx),2);
	b(2*id(~fdx))   = +c(id(~fdx),3)./c(id(~fdx),1);


	A = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),2*nxc,2*nxc);

	% boundary condition at left
	% only one boundary has to be specified for first order odes
	%[v, p, ~, set] = bcfun(xi(1),[],ccdx);
	[v, p, ~, set] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});

	nb = 0;
	if (set) %~isempty(v))
		% p*f(0) + (1-p)*f(0)'' = v
		row = 2*nxc; % -1
		if (0 ~= c(1,2))
			A(row,1)  = p*obj.exp(-0.5*r(1)*dx(1)) + (1-p)*r(1)*obj.exp(-0.5*r(1)*dx(1));
			A(row,2) = p;
			b(row)   = v;
		else
			nbuf = nbuf+1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 1;
			Abuf(nbuf,3) = p*(-0.5*dx(1)) + (1-p);
			nbuf = nbuf+1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 2;
			Abuf(nbuf,3) = p;
			A(row,1) = p*(-0.5*dx(1)) + (1-p);
			A(row,2) = p;
			b(row)   = v;
		end
		nb           = 1;
	end

	%[v, p, ~, set] = bcfun(xi(2),[],ccdx);
	[v, p, ~, set] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});

	if (set) %~isempty(v))
		row=2*nxc; % -1
		if (0~=c(end,2))
			A(row,2*nxc-1) = (       p*obj.exp(+0.5*r(nxc)*dx(nxc)) ...
					     + (1-p)*r(nxc)*obj.exp(+0.5*r(nxc)*dx(nxc)) ...
					   );
			A(row,2*nxc) = p;
			b(row)       = v;
		else
			nbuf         = nbuf + 1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 2*nxc-1;
			Abuf(nbuf,3) = ( p*0.5*dx(nxc) + (1-p) );
			nbuf         = nbuf + 1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 2*nxc;
			Abuf(nbuf,3) = p;
			A(row,2*nxc-1) = ( p*0.5*dx(nxc) + (1-p) );
			A(row,2*nxc)   = p;
			b(row)       = v;
		end
		nb     = nb+1;
	end
	A_ = sparse(Abuf(:,1),Abuf(:,2),Abuf(:,3),2*nxc,2*nxc);
	Abuf(:,1:2) = Abuf(:,1:2) + npi(edx)-1 + obj.npii(1,cdx)-1;
	obj.Abuf   = [obj.Abuf; Abuf];
	obj.b(obj.npii(edx,cdx):obj.npii(edx+1,cdx)-1,1)   = b;

	if (1 ~= nb)
		error('need boundary condition at exactly one end');
	end
end % assemble1_A

