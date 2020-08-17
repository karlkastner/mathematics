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
% TODO use buffers instead of sparse
function [out] = bvp2c(odefun,bcfun,ifun,xi,nx,varargin)

	opt = bvp2_check_arguments(varargin{:});

	% number of edges in graph (channels in network)
	nc = length(odefun);

	% solution for edges in graph (channel in network)
	out = struct();

	% grid for each each edge (channel in network)
	for cdx=1:nc 
		% segment end points
		out(cdx).x = mesh1(xi(cdx,:),nx(cdx),opt.xs);
	
		% segment mid points
		out(cdx).xc = mid(out(cdx).x);
	
		% segment lengths
		out(cdx).dx  = diff(out(cdx).x);
	end

	% number of segments
	nxc = nx-1;
	
	% number of parallel-coupled odes, identical for all edges (channels)
	oo  = feval(odefun{1});

	% number of equations per segment
	neq = length(oo);

	% start index of segment end-points
	ni    = zeros(neq+1,nc);
	ni(1,:) = 1;
	% start index of segment mid-points
	nci    = zeros(neq+1,nc);
	nci(1,:) = 1;
	% start index of segment-mid points, separated parts of ode solution
	npi    = zeros(neq+1,nc);
	npi(1,:) = 1;

	% number of equations per ode
	for cdx=1:neq
		% oder of ode
		switch (oo(cdx))
		case {-1} % special for Q0 coupling
			ni(cdx+1,:)  = ni(cdx,:)  + 1;
			nci(cdx+1,:) = nci(cdx,:) + 1;
			npi(cdx+1,:) = npi(cdx,:) + 1;
		case {1} % first order
			% homogeneous and inhomogeneous part
			ni(cdx+1,:)  = ni(cdx,:)   + nxc+1;
			nci(cdx+1,:)  = nci(cdx,:) + nxc;
			npi(cdx+1,:) = npi(cdx,:)  + 2*nxc;
		case {2} % second order
			% homogeneous left-going, right-going and inhomogeneous part
			ni(cdx+1,:)  = ni(cdx,:)  + nxc+1;
			nci(cdx+1,:) = nci(cdx,:) + nxc;
			npi(cdx+1,:) = npi(cdx,:) + 3*nxc;
		otherwise
			error('');
		end
	end % for neq

	% indices into global discretization matrix
	npii   = npi + [0,npi(end,1:end-1)-1];

	% initial value of ypm
	% complex amplitude of the left and right going wave at segment mid points
	% TODO should reassembly not take place inside iteration?
	ypm    = zeros(npii(end,end)-1,1);
	if (~isempty(ifun))
	    for cdx=1:nc
		yi      = feval(ifun{cdx},out(cdx).x);
		for edx = 1:neq
			% assign initial values into the inhomogeneous part
			% this does not matter, as yi is reassembled later
			% the inhomogeneous part is always the second part
			yi_ = yi(ni(edx,cdx):ni(edx+1,cdx)-1);
			switch (oo(edx))
			case {-1}
				ypm(npi(edx),cdx) = yi_;
			case {1}
				yci = mid(yi_);
				ypm(npi(edx,cdx):2:npi(edx+1,cdx)-2) = yci;
			case {2}
				yci = mid(yi_);
				ypm(npi(edx,cdx):3:npi(edx+1,cdx)-2) = yci;
			end % switch
		end % for edx
	    end % for cdx
	end % ~isempty(ifun)

	% solve non-linear system by picard iteration
	[ypm, out.cflag, out.kiter] = picard(@bvp2c_solve,ypm,opt.sopt);

	% interpolate solution to segment end points (grid points)
	inner2outer_();

function ypm = bvp2c_solve(ypm)
	[AA,bb,out] = bvp2c_assemble(out,ypm,odefun,bcfun,xi,neq,nxc,oo,ni,nci,npi,npii);

	% balance
	if (isfield(opt,'balance') && opt.balance)
		s  = 1./sqrt(abs(diag(AA)));
		AA = diag(s)*AA.*diag(s);
		bb = s.*bb;
	end

	% solve
	ypm    = (AA \ bb);
end % bvp2c_solve

	function inner2outer_()
	    for cdx=1:nc
		out(cdx).y = zeros(ni(end,cdx)-1,1);
		for edx=1:neq
			r = out(cdx).ll(:,1,edx);
			switch (oo(edx))
			case {-1}
				y_ = ypm(npii(edx,cdx));
			case {1}
			    if (0 ~= r(1,1,edx))
				y_ = (   ypm(npii(edx,cdx)  :2:npii(edx+1,cdx)-2) ...
				       + ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1) ...
				     );
			    else
				% degenerated, linear function
				y_ = ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1); % or +0, -2?
			    end
				y_ = inner2outer(y_);
			case {2}
				% TODO use bvp1c/2c to expand with exact exponential
				y_ = (   ypm(npii(edx,cdx)  :3:npii(edx+1,cdx)-3) ...
				       + ypm(npii(edx,cdx)+1:3:npii(edx+1,cdx)-2) ...
				       + ypm(npii(edx,cdx)+2:3:npii(edx+1,cdx)-1) ...
			             );
				y_ = inner2outer(y_);
		        otherwise
				error('here');
		     	end % switch mm
			out(cdx).y(ni(edx,cdx):ni(edx+1,cdx)-1) = y_;
		end % for edx
	    end % for cdx
	end % inner2outer_

end % function bvp2c

