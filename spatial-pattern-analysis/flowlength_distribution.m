% Mon 25 Mar 15:57:59 CET 2024
% Karl Kastner, Berlin
%
% for stationary pattern dh/dt = 0 and diffusive flow
% with source (precipitation) and sink (infiltration)
%
% 0 = I + R + eh*(d2h/dx2 + d2h/dy2)
%
%
%
% uses manhatten distance
%
function [P, Pij] = flowlength_distribution(h,I,R,ew,dmax)
    % local 3x3 4-neighbour indices (south,west,north,east)
    ni = [-1, 0, 0, 1];
    nj = [0,-1, 1, 0];

    n   = size(h);

	
    mean(I,'all')+mean(R,'all')
pause

    visited = zeros(n,'logical');
    Pij = zeros(dmax,n(1),n(2));
    
    % for each cell
    for idx=1:n(1)
	[idx]
        for jdx=1:n(2)
    		flow_distribution_(idx,jdx)
        end
    end
    % the mean distribution is just the average of all cells
    P = mean(mean(Pij,3),2);

    function flow_distribution_(idx,jdx)
	if (~visited(idx,jdx))
		% has to be set first to prevend infinite recursion
		visited(idx,jdx) = true;
		% for each neighbour
		% get flow distribution of upstream cell
		for ndx=1:length(ni)
			% cicular boundary conditions
			in = mod(idx+ni(ndx)-1,n(1))+1;
			jn = mod(jdx+nj(ndx)-1,n(2))+1;
			% if neigbour is an upstream cell 
			if (h(in,jn) > h(idx,jdx))
				flow_distribution_(in,jn);
				% increase flowlength by one for water coming from neighbouring cell
				% weight = flow from neighbour to this cell / total water height in neibouring cell
				w = ew*(h(in,jn)-h(idx,jdx))/h(in,jn);
				Pij(2:end,idx,jdx) = Pij(2:end,idx,jdx) + w*Pij(1:end-1,in,jn);
			end
		end % for ndx 
		% precipitation into current cell has flowlength 0
    		if (isscalar(R))
			Pij(1,idx,jdx) = R;
		else
			Pij(1,idx,jdx) = R(idx,jdx);
		end
		% subtract amount of water that is infiltrating,
		% water volumes of each flowlength is reduced proportionally
		w = 1/sum(Pij(:,idx,jdx));
		% note : the limitation seems to be necessary to account for the nonzero dh/dt
		Pij(:,idx,jdx) = max(0,1+I(idx,jdx)*w)*Pij(:,idx,jdx);
		% TODO there seems to be an issue here, without flow, everything infiltrates with length 0 and all elements in Pij = 0
        end % if visited
    end % flow_distribution_

end % flow_distribution

