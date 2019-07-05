% 2016-03-15 12:26:55.438222422 +0100
function tV = interp_sn_(idx,sdx,sS,sN,tS,sV,tN,tV,ns,order,L_max2)
		if (~isempty(sdx))
			% get 2-point segments from point track
			seg = [cvec(sdx), min(cvec(sdx+1),ns)];
			% discard segments that are too long
			L2   = (sS(seg(:,2))-sS(seg(:,1))).^2 + (sN(seg(:,2))-sN(seg(:,1))).^2;
			sdx_ = L2 <= L_max2;
			seg  = seg(sdx_,:);

			% select only segments that are convex in this point
			% TODO, allow for extrapoltaion
			c    = (sN(seg(:,2)) - tN(idx))./(sN(seg(:,2))-sN(seg(:,1)));
			sdx_ = (c >= 0) & (c <= 1);
			seg  = seg(sdx_,:);
			if (size(seg,1) > 0)
				c    = c(sdx_);
		%		sdx  = sdx(sdx_); 
				% interpolate points to the N-coordinate
				S_ = c.*sS(seg(:,1)) + (1-c).*sS(seg(:,2));
		
				% determine closest convex segment ahead of s0
				switch (order)
				case {-1} % mark only for test
					tV(idx,:) = length(S_);
				case {0} % nearest neighbour
					% this is linear in N and nearest neighbour along S
					[S_min mdx] = min(abs(S_-tS(idx)));
					tV(idx,:) = c(mdx)*sV(seg(mdx,1),:) + (1-c(mdx))*sV(seg(mdx,2));
				case {1} % linear
					dS  = S_ - tS(idx);
					ldx = find(dS <= 0);
					rdx = find(dS >= 0);
					% fall back to constant interpolation if either ldx or sdx is empty
					if (isempty(ldx) || isempty(rdx))
						[S_min mdx] = min(abs(S_-tS(idx)));
						tV(idx,:) = c(mdx)*sV(seg(mdx,1),:) + (1-c(mdx))*sV(seg(mdx,2));
					else
					[dSl mdx] = min(abs(dS(ldx)));
					ldx = ldx(mdx);
					Vl  = c(ldx)*sV(seg(ldx,1),:) + (1-c(ldx))*sV(seg(ldx,2),:);
					[dSr mdx] = min(abs(dS(rdx)));
					rdx = rdx(mdx);
					Vr  = c(rdx)*sV(seg(rdx,1),:) + (1-c(rdx))*sV(seg(rdx,2),:);
					c       = (S_(rdx) - tS(idx))/(S_(rdx) - S_(ldx));
					tV(idx,:) = c*Vl + (1-c)*Vr;
					end
				otherwise
					error('not yet implemented');
				end
			end % ~isempty(seg)
		end % ~isempty(sdx)

