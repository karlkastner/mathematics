% Di 10. sNov 17:16:14 CET 2015
% Karl Kastner, Berlin
%
%% interpolate along streamwise coordinates
%% This gives similar result to setting aspect ratio for sN to infinity,
%% but not quite,as the input point set is not dense (scale for sN to infinity does not work)
% TODO at bifurcations: follow deeper channel, or follow channel that minimises difference in cross section
% TODO special treatment to outside points
% add artifical lines from nmax/w to 0.5 and nmin/w to -0.5

function tV = interp_sn(sS,sN,sV,tS,tN,sS_range,L_max,order,X,Y)

	% create search object along S
	%sobj = KDTreesSearcher(sS);
	sobj = ExhaustivesSearcher(sS);

	% allocate memory
	ns      = length(sS);
	nt      = length(tS);
	tV     = NaN(n,size(sV,2));
	L_max2 = L_max*L_max;
	
	% for each target point
	for idx=1:nt
		% get all points within range
		% TODO, only in the same segment
		sdx = rangesearch(sobj,tS(idx,:),sS_range);
		sdx = (1:nt)';
		tV = interp_sn_(idx,sdx,sS,SN,tS,sV,tN,tV,ns,order,L_max2);
%	if (~isempty(sdx))
%		% get 2-point segments from point track
%		seg = [cvec(sdx), min(cvec(sdx+1),ns)];
%		% discard segments that are too long
%		L2 = (sS(seg(:,2))-sS(seg(:,1))).^2 + (sN(seg(:,2))-sN(seg(:,1))).^2;
%		sdx_ = L2 <= L_max2;
%		seg  = seg(sdx_,:);
%
%		% select only segments that are convex in this point
%			% TODO, allow for extrapoltaion
%		c    = (sN(seg(:,2)) - tN(idx))./(sN(seg(:,2))-sN(seg(:,1)));
%		sdx_ = (c >= 0) & (c <= 1);
%		seg  = seg(sdx_,:);
%		if (size(seg,1) > 0)
%		c    = c(sdx_);
%%		sdx  = sdx(sdx_); 
%		% interpolate points to the N-coordinate
%		S_ = c.*sS(seg(:,1)) + (1-c).*sS(seg(:,2));
%
%		% determine closest convex segment ahead of s0
%		switch (order)
%				case {-1} % mark only for test
%					tV(idx,:) = length(S_);
%		case {0} % nearest neighbour
%			% this is linear in N and nearest neighbour along S
%			[S_min mdx] = min(abs(S_-tS(idx)));
%			tV(idx,:) = c(mdx)*sV(seg(mdx,1),:) + (1-c(mdx))*sV(seg(mdx,2));
%		case {1} % linear
%			dS  = S_ - tS(idx);
%			ldx = find(dS <= 0);
%			rdx = find(dS >= 0);
%			% fall back to constant interpolation if either ldx or sdx is empty
%					if (isempty(ldx) || isempty(rdx))
%						[S_min mdx] = min(abs(S_-tS(idx)));
%						tV(idx,:) = c(mdx)*sV(seg(mdx,1),:) + (1-c(mdx))*sV(seg(mdx,2));
%					else
%			[dSl mdx] = min(abs(dS(ldx)));
%			ldx = ldx(mdx);
%			Vl  = c(ldx)*sV(seg(ldx,1),:) + (1-c(ldx))*sV(seg(ldx,2),:);
%			[dSr mdx] = min(abs(dS(rdx)));
%			rdx = rdx(mdx);
%			Vr  = c(rdx)*sV(seg(rdx,1),:) + (1-c(rdx))*sV(seg(rdx,2),:);
%			c       = (S_(rdx) - tS(idx))/(S_(rdx) - S_(ldx));
%			tV(idx,:) = c*Vl + (1-c)*Vr;
%			end
%		otherwise
%			error('not yet implemented');
%		end
%		end % ~isempty(seg)
%		end % ~isempty(sdx)
	end % for idx

end % function

