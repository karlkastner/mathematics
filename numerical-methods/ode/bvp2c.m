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
% TODO make classe
% TODO allow for resolve (do not rerun through init process)
function [out] = bvp2c(odefun,bcfun,ifun,xi,nx,junction_condition,dischargeisvariable,varargin)

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
	ni       = zeros(neq+1,nc);
	ni(1,:)  = 1;
	% start index of segment mid-points
	nci      = zeros(neq+1,nc);
	nci(1,:) = 1;
	% start index of segment-mid points, separated parts of ode solution
	npi      = zeros(neq+1,nc);
	npi(1,:) = 1;

	% number of equations per ode
	for cdx=1:neq
		% oder of ode
		switch (oo(cdx))
		case {1} % first order
			% homogeneous and inhomogeneous part
			ni(cdx+1,:)  = ni(cdx,:)   + rvec(nxc)+1 + dischargeisvariable;
			nci(cdx+1,:) = nci(cdx,:)  + rvec(nxc)   + dischargeisvariable;
			npi(cdx+1,:) = npi(cdx,:)  + 2*rvec(nxc) + dischargeisvariable;
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

	% initial value of ypm
	% complex amplitude of the left and right going wave at segment mid points
	ypm    = zeros(npii(end,end)-1,1);
	if (~isempty(ifun))
	    for cdx=1:nc
		yi      = feval(ifun{cdx},out(cdx).x);
		if (length(yi) == npi(end,cdx)-1)
				ypm(npii(1,cdx):npii(end,cdx)-1) = yi;
		else
		for edx = 1:neq
			% assign initial values into the inhomogeneous part
			% this does not matter, as yi is reassembled later
			% the inhomogeneous part is always the second part
			if (length(yi) == nci(end,cdx)-1)
				yi_ = yi(nci(edx,cdx):nci(edx+1,cdx)-1);
			else
				yi_ = yi(ni(edx,cdx):ni(edx+1,cdx)-1);
			end
			switch (oo(edx))
			case {1}
				if (~dischargeisvariable)
					if (length(yi) ~= nci(end,cdx)-1)
						yi_ = mid(yi_);
					end
					ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1) = yi_;
				else
					% discharge
					ypm(npii(edx+1,cdx)-1) = yi_(end);

					% level
					if (length(yi) ~= nci(end,cdx)-1)
						yi_ = mid(yi_(1:end-1));
					end
					ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-2) = yi_;
				end
			case {2}
				if (length(yi) ~= nci(end,cdx)-1)
					yi_ = mid(yi_);
				end
				ypm(npii(edx,cdx)+1:3:npii(edx+1,cdx)-1) = yi_;
			end % switch
		   end % if
		end % for edx
	    end % for cdx
	end % ~isempty(ifun)
%figure(1)
%clf
%	plot(abs(ypm))
%plot([real(ypm),imag(ypm)])

	% solve non-linear system by picard iteration
	[ypm, out(1).cflag, out(1).kiter] = picard(@bvp2c_solve,ypm,opt.sopt);

	% interpolate solution to segment end points (grid points)
	reconstruct();

function ypm = bvp2c_solve(ypm)
	[AA,bb,out] = bvp2c_assemble(out,ypm,odefun,bcfun,xi,neq,nxc,oo,ni, ...
			  nci,npi,npii,junction_condition,dischargeisvariable,opt.bcarg);

	% balance
	if (isfield(opt,'balance') && opt.balance)
		s  = 1./sqrt(abs(diag(AA)));
		AA = diag(s)*AA.*diag(s);
		bb = s.*bb;
	end

	% solve
	ypm    = (AA \ bb);
%	ypm = pinv(full(AA))*bb;
end % bvp2c_solve

	% TODO reconstruct with exp instead of using inner2outer
	function reconstruct()
	    for cdx=1:nc
	        yc = zeros(nci(end,cdx)-1,1);
		if (opt.reconstruct_y)
			y = zeros(ni(end,cdx)-1,1);
		end
		for edx=1:neq
			c = out(cdx).cc;
			switch (oo(edx))
			case {1}
			    if (~dischargeisvariable)
				% TODO, check for each row
			     if (0 ~= c(1,2,edx))
				% this does not happen for the
				yc_ = (   ypm(npii(edx,cdx)  :2:npii(edx+1,cdx)-2) ...
				       + ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1) ...
				     );
			     else
				% degenerated, linear function
				yc_ = ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-1); % or +0, -2?
			     end
				if (opt.reconstruct_y)
					y_ = inner2outer(yc_);
				end
			    else
			     if (0 ~= c(1,2,edx))
				yc_ = (  ypm(npii(edx,cdx)  :2:npii(edx+1,cdx)-3) ...
				       + ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-2) ...
				      );
			     else
				% degenerated, linear function
				yc_ = ypm(npii(edx,cdx)+1:2:npii(edx+1,cdx)-2); % or +0, -2?
			     end
				if (opt.reconstruct_y)
					y_  = [inner2outer(yc_); ypm(npii(edx+1,cdx)-1)];
				end
				yc_ = [yc_; ypm(npii(edx+1,cdx)-1)];
			    end
			case {2}
				% TODO use bvp1c/2c to expand with exact exponential
				yc_ = (  ypm(npii(edx,cdx)  :3:npii(edx+1,cdx)-3) ...
				       + ypm(npii(edx,cdx)+1:3:npii(edx+1,cdx)-2) ...
				       + ypm(npii(edx,cdx)+2:3:npii(edx+1,cdx)-1) ...
			             );
				if (opt.reconstruct_y)
					y_ = inner2outer(yc_);
				end
		        otherwise
				error('here');
		     	end % switch mm
			yc(nci(edx,cdx):nci(edx+1,cdx)-1) = yc_;
			if (opt.reconstruct_y)
				y(ni(edx,cdx):ni(edx+1,cdx)-1) = y_;
			end
		end % for edx
	    	out(cdx).yc   = yc;
		out(cdx).y    = y;
		out(cdx).ypm  = ypm(npii(1,cdx):npii(end,cdx)-1);
	    end % for cdx
	end % reconstruct

end % function bvp2c

