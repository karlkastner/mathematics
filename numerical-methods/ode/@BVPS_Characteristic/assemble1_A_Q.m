% Sat 18 Apr 22:02:13 +08 2020
function assemble1_A_Q(obj,cdx,edx)
	xi     = obj.xi(cdx,:);
	xc     = obj.out(cdx).xc;
	dx     = obj.out(cdx).dx;
	nxc    = obj.nxc(cdx);

	% only use last 3 coefficients
	%c  = cc(:,end-2:end,ccdx); 

	% c1 z' + c2 z + c3 Q0 + c4 == 0
	c  = obj.out(cdx).cc(:,:,edx); 
	% TODO explot that c==0 when r==0
	fdx = (c(:,2) ~= 0);

	% root for homogeneous part
	r  = obj.out(cdx).ll(:,1,edx);

	% match at rightend at segment i with right end at segment i+1
	nxc = size(c,1);

	Abuf = zeros(6*(nxc-1)+4,3);
	b    = zeros(2*nxc+1,1);

	% first row : left boundary condition
	% row 2:2:2*nxc : ode
	% row 3:2:2*nxc : matching condition (continuity)

	nbuf = 0;
	% continuities at segment interfaces
	id = (1:nxc-1)';
	Abuf(nbuf+id,1) = 2*id+1;
	Abuf(nbuf+id,2) = 2*id-1; % 1,3,5, ...
	% match values at right and left end of segment j and j+1
	% degenerated, without damping,
	% the piecewise solution is a linear function
	Abuf(nbuf+id,3)   = (   fdx(1:end-1).*( -obj.exp(+0.5*r(id).*dx(id))) ...
			      + (1-fdx(1:end-1)).*-0.5.*dx(id));
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+id,1) = 2*id+1;
	Abuf(nbuf+id,2) = 2*id;
	Abuf(nbuf+id,3) = -1;
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+id,1) = 2*id+1;
	Abuf(nbuf+id,2) = 2*id+1;
	Abuf(nbuf+id,3)  = ( fdx(2:end).*+obj.exp(-0.5*r(id+1).*dx(id+1)) ...
			         + (1-fdx(2:end)).*-0.5.*dx(id+1) );
	nbuf = nbuf+nxc-1;
	Abuf(nbuf+id,1) = 2*id+1;
	Abuf(nbuf+id,2) = 2*id+2;
	Abuf(nbuf+id,3) = +1;
	nbuf = nbuf+nxc-1;

	% ode, to determine the inhomogeneous part and Q0
	id              = (1:nxc)';
	Abuf(nbuf+id,1) = 2*id;
	Abuf(nbuf+id,2) = 2*id + (fdx-1); % 1,3,5
	Abuf(nbuf+id,3) = c(:,1);
	nbuf = nbuf+nxc;

	% TODO c2 z (not necessary here, but for completeness)
	Abuf(nbuf+id,1) = 2*id;
	Abuf(nbuf+id,2) = 2*nxc+1;
	Abuf(nbuf+id,3) = c(:,3);
	nbuf = nbuf+nxc;
	b(2*id) = -c(:,4);

	% boundary condition at left end
	% only one boundary has to be specified for first order odes
	%[v, p, ~, type] = bcfun(xi(1),[],ccdx);
	[v, p, ~, type] = obj.bcfun(cdx,1,edx,obj.opt.bcarg{:});
	row = 1;
	switch (type)
	case {'z'}
		% p*f(0) + (1-p)*f(0)'' = v
		if (0~=c(1,2))
			A(row,1) = p(1)*obj.exp(-0.5*r(1)*dx(1)) + (1-p)*r(1)*obj.exp(-0.5*r(1)*dx(1));
			A(row,2) = p;
			b(row)   = v;
		else
			nbuf        = nbuf+1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 1;
			Abuf(nbuf,3) = p(1)*(-0.5*dx(1)) + (1-p(1));
			nbuf = nbuf+1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 2;
			Abuf(nbuf,3) = p(1);
			b(row)       = v;
		end
		nQ = 0;
	case {'Q'}
		%row  = 2*nxc+1;
		nbuf = nbuf+1;
		Abuf(nbuf,1) = row;
		Abuf(nbuf,2) = 2*nxc+1;
		Abuf(nbuf,3) = 1;
		nbuf         = nbuf+1;
		Abuf(nbuf,1) = row;
		Abuf(nbuf,2) = 2*nxc+1;
		Abuf(nbuf,3) = 0; % dummy, for non-zero index in Abuf
		b(row) = v; 
		nQ = 1;
	case {''}
		nQ = 0;
		% nothing to do
	otherwise
		error('here');
	end

	%[v, p, ~, type] = bcfun(xi(2),[],ccdx);
	[v, p, ~, type] = obj.bcfun(cdx,2,edx,obj.opt.bcarg{:});
	row=2*nxc+1;
	switch (type)
	case {'z'}
		if (0~=c(end,2))
			A(row,2*nxc-1) = (       p(1)*obj.exp(+0.5*r(nxc)*dx(nxc)) ...
					     + (1-p(1))*r(nxc)*obj.exp(+0.5*r(nxc)*dx(nxc)) ...
					 );
			A(row,2*nxc)   = p(1);
			b(row)         = v;
		else
			nbuf         = nbuf + 1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 2*nxc-1;
			Abuf(nbuf,3) = ( p(1)*0.5*dx(nxc) + (1-p(1)) );
			nbuf         = nbuf + 1;
			Abuf(nbuf,1) = row;
			Abuf(nbuf,2) = 2*nxc;
			Abuf(nbuf,3) = p;

			b(row)       = v;
		end
	case {'Q'}
		if (1 == nQ)
			error('need, discharge can only be specified once, need at least level at one end');
		end
		nbuf = nbuf+1;
		Abuf(nbuf,1) = row;
		Abuf(nbuf,2) = 2*nxc+1;
		Abuf(nbuf,3) = 1;
		nbuf = nbuf+1;
		Abuf(nbuf,1) = row;
		Abuf(nbuf,2) = 2*nxc+1;
		Abuf(nbuf,3) = 0; % dummy
		b(row) = v; 
	case {''}
		% nothing to do
	otherwise
		error('here');
	end

	Abuf(:,1:2) = Abuf(:,1:2) + obj.npi(edx,cdx)-1 + obj.npii(1,cdx)-1;
	obj.Abuf   = [obj.Abuf; Abuf];
	obj.b(obj.npii(edx,cdx):obj.npii(edx+1,cdx)-1,1)   = b;
end % assemble1_A_Q

