% Fr 15. Jan 09:55:03 CET 2016
% Karl Kastner, Berlin
%
%% rotation matrix from angle
function R = rot2(alpha,idx,jdx,n)
	if (nargin()<2)
		idx=1;
	end
	if (nargin()<3)
		jdx=2;
	end
	if (nargin() < 3)
		n=max(idx,jdx);
	end
%	idx = min(idx_,jdx_);
%	jdx = max(idx_,jdx_);
	R = speye(n);
	c = cos(alpha);
	s = sin(alpha);
	R(idx,idx) =  c;
	R(jdx,jdx) =  c;
	R(idx,jdx) = -s; % signs were swapped
	R(jdx,idx) = +s;
end

