% Wed 26 Aug 14:58:23 +08 2020

function init(obj)
	obj.check_arguments();

	% number of edges in graph (channels in network)
	obj.nc = size(obj.xi,1);

	% solution for edges in graph (channel in network)
	obj.out = struct();

	% grid for each each edge (channel in network)
	for cdx=1:obj.nc 
		% segment end points
		obj.out(cdx).x = mesh1(obj.xi(cdx,:),obj.nx(cdx),obj.opt.xs);
	
		% segment mid points
		obj.out(cdx).xc = mid(obj.out(cdx).x);
	
		% segment lengths
		obj.out(cdx).dx  = diff(obj.out(cdx).x);
	end

	% number of segments in each channel
	nxc = obj.nx-1;
	
	% number of parallel-coupled odes, identical for all edges (channels)
	oo  = feval(obj.odefun,1);

	% number of equations per segment
	neq = length(oo);

	% start index of segment end-points
	ni       = zeros(neq+1,obj.nc);
	ni(1,:)  = 1;
	% start index of segment mid-points
	nci      = zeros(neq+1,obj.nc);
	nci(1,:) = 1;
	% start index of segment-mid points, separated parts of ode solution
	npi      = zeros(neq+1,obj.nc);
	npi(1,:) = 1;

	% number of equations per ode
	for cdx=1:neq
		% oder of ode
		switch (oo(cdx))
		case {1} % first order
			% homogeneous and inhomogeneous part
			ni(cdx+1,:)  = ni(cdx,:)   + rvec(nxc)+1 + obj.opt.dischargeisvariable;
			nci(cdx+1,:) = nci(cdx,:)  + rvec(nxc)   + obj.opt.dischargeisvariable;
			npi(cdx+1,:) = npi(cdx,:)  + 2*rvec(nxc) + obj.opt.dischargeisvariable;
		case {2} % second order
			% homogeneous left-going, right-going and inhomogeneous part
			ni(cdx+1,:)  = ni(cdx,:)  + rvec(nxc)+1;
			nci(cdx+1,:) = nci(cdx,:) + rvec(nxc);
			npi(cdx+1,:) = npi(cdx,:) + 3*rvec(nxc);
		otherwise
			error('');
		end
	end % for neq

	% indices into global discretization matrix
	npii   = npi + [0,cumsum(npi(end,1:end-1)-1)];

	obj.neq = neq;
	obj.nxc = nxc;
	obj.oo  = oo;
	obj.ni  = ni;
	obj.nci = nci;
	obj.npi = npi;
	obj.npii = npii;
end % init

