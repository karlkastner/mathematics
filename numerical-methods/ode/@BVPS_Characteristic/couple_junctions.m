% Sat 28 Oct 14:43:06 CEST 2017
% TODO, treatment of degenerated first order solution at junctions
%
function [AAA,bb] = couple_junctions(obj,AAA,bb)
    if (~obj.opt.dischargeisvariable)
		error('here')
    end
	xi = obj.xi;
	npii = obj.npii;

    % mid point of channels (for determining direction)
    ximid = 0.5*(xi(:,1)+xi(:,2));

    % TODO, make general, this is semi-customized for river tide computation
    % for each junction
    for jdx=1:length(obj.jfun)

	% channel direction
	% dir = sign(xi(:,2)-xi(:,1));
    
    	% cid : channel indices
	% eid : endpoint indices {1,2},
	% p   : scales, for matching level of the frequency compontents,
	%               these are 1/(iow) 
    	[cid, eid, p] = feval(obj.jfun{jdx});

	% fetch end point properties
	% segment length
	dx = zeros(length(cid),1);

	% eigenvalues
	l  = zeros(length(cid),2,obj.neq);
	% indices
	rid_ = zeros(obj.neq,length(cid));
	col_ = zeros(obj.neq,length(cid));
	for cdx=1:length(cid)
		switch (eid(cdx))
		case {1}
			dx(cdx)     = obj.out(cid(cdx)).dx(1);
			l(cdx,:,:)  = obj.out(cid(cdx)).ll(1,:,:);
			rid_(:,cdx) = npii(1:end-1,cid(cdx));
			col_(:,cdx) = npii(1:end-1,cid(cdx));
		case {2}
			dx(cdx)     = obj.out(cid(cdx)).dx(end);
			l(cdx,:,:)  = obj.out(cid(cdx)).ll(end,:,:);
			rid_(:,cdx) = npii(2:end,cid(cdx))-1;  
			% TODO this is a quick and dirty hack
			% nb : tidal shift by 3, because there are 3 parts
			%      mean shift by 3, because there are 2 parts and the mean discharge
			col_(:,cdx) = npii(2:end,cid(cdx)) - 3.*(1==cvec(obj.oo)) - 3.*(2==cvec(obj.oo));
			%rid_(:,cdx) = rid_(:,cdx) - 1;
		otherwise
			error('here');
		end
    		dir(cdx) = sign(xi(cdx,eid(cdx)) - ximid(cid(cdx)));
	end % for cdx

    	% for each parallel equation (frequency component)
    	for edx=1:obj.neq
    		% row indices
		rid = rid_(edx,:)';
		col = col_(edx,:)';

		switch (obj.oo(edx))
		case {1} % match mean frequency component
			% condition for first connecting channel :
    			% conservation of (tidally averaged) discharge
			% sum dir*Q0_i  == 0 (1 equation)
			% reset row
			AAA(rid(1),:) = 0;
			bbb(rid(1))   = 0;
			for cdx=1:length(cid)
				% mean discharge Q0 has only 1 variable
				%AAA(rid(1),rid(cdx)) = dir(cdx);
				% here, the column indices are for the tidally averaged discharge,
				% not of the tidally averaged water level (rid)
				% TODO avoid magic number for index "2"
				AAA(rid(1),npii(2,cid(cdx))-1) = dir(cdx);
			end % for cdx
				
			% condition for remaining channels meeting at junction:
			% continuty of surface elevation (equal tidally averaged water level)
    			% z0_i - z0_1 == 0, i>1 (n-1 equations)
			% the water level does not depend on channel direction,
			% so dir is not factored in
    			for cdx=2:length(cid)
				% reset row
				AAA(rid(cdx),:) = 0;
				bbb(rid(cdx))   = 0;
			% note : solution for water level is always degenerated
				% homogeneous part	
% TODO check cc (!)
				if (0) %abs(l(1,1,edx)) ~= 0)
					% this actually never happens
					AAA(rid(cdx),col(1))     = -obj.exp(dir(1)*0.5*l(1,1,edx).*dx(1));
				else
					AAA(rid(cdx),col(1))     = -dir(1)*0.5*dx(1);
				end
				% inhomogeneous part
				AAA(rid(cdx),col(1)+1)   = -1;

				% homogeneous part	
% TODO check cc (!)
				if (0) %abs(l(cdx,1,edx)) ~= 0)
					% this actually never happens
					AAA(rid(cdx),col(cdx))   = +obj.exp(dir(cdx)*0.5*l(cdx,1,edx).*dx(cdx));
				else
					AAA(rid(cdx),col(cdx))   = +dir(cdx)*0.5*dx(cdx);
				end
				% inhomogeneous part
				AAA(rid(cdx),col(cdx)+1) = +1;
    			end % for cdx
		case {-1}
			% TODO deprecated
			% nothing to do (Q0 is matched for 1==edx)
		case {2}
			% frequency components (tide)

			% condition for first connecting channel :
			% conservation of tidal discharge of e-th frequency component
			% sum s Qt_i == 0, i = 1 .. n

			% reset row
    			AAA(rid(1),:) = 0;
			% no external inflow/outflow at junction
    			bbb(rid(1))   = 0;
    
    			for cdx=1:length(rid)
    				% discharge consits of 3 unknowns (columns)
				% left going part (Q-)
    				AAA(rid(1),col(cdx))   = dir(cdx)*obj.exp(dir(cdx)*0.5*l(cdx,1,edx)*dx(cdx));
				% inhomogeneous part
    				AAA(rid(1),col(cdx)+1) = dir(cdx);
				% right going part (Q+)
				AAA(rid(1),col(cdx)+2) = dir(cdx)*obj.exp(dir(cdx)*0.5*l(cdx,2,edx)*dx(cdx));
    			end % for cdx

			% condition for remaining connecting channels :
			% equal (oscillation) of water level of e-th frequency component
			%    zt_i - zt_1 == 0
			% => 1/(iow_i) dQt_i/dx - 1/(iow_1) dQ_1/dx == 0
			% => 1/(iow) (l_i^- Q_i^- + l_i^+ Q_i^+) - 1/(iow_1) (l_i^- Q_i^- + l_i^+ Q_i^+) == 0
			% p = 1/iow, 1/io can be factored out, but 1/w not
			% sign of channel direction is not factored in, as the water level is not a directional vector quantity
			for cdx=2:length(rid)
				% reset row
    				AAA(rid(cdx),:) = 0;
    				bbb(rid(cdx))   = 0;

				% left going
    				AAA(rid(cdx),col(1))     = -p(1)*l(1,1,edx)*obj.exp(dir(1)*0.5*l(1,1,edx)*dx(1));
				% derivative of inhomogeneous part is zero
				% right going
				AAA(rid(cdx),col(1)+2)   = -p(1)*l(1,2,edx)*obj.exp(dir(1)*0.5*l(1,2,edx)*dx(1));

				% left going
    				AAA(rid(cdx),col(cdx))   = +p(cdx)*l(cdx,1,edx)*obj.exp(dir(cdx)*0.5*l(cdx,1,edx)*dx(cdx));
				% derivative of inhomogeneous part is zero
				% right going
				AAA(rid(cdx),col(cdx)+2) = +p(cdx)*l(cdx,2,edx)*obj.exp(dir(cdx)*0.5*l(cdx,2,edx)*dx(cdx));
			end % for cdx
		otherwise
			error('here');
		end % switch oo(edx)
	end % for edx (each frequ ency component)
    end % for cdx (each junction)
end % function coupling_condition

