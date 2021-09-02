% Fr 15. Jan 09:55:03 CET 2016
% Karl Kastner, Berlin
%
%% rotation matrix from angle
% function R = rot2(alpha,idx,jdx,n)
function R = rot2(alpha_rad,idx,jdx,n)
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
if (~issym(alpha_rad))
	R = speye(n);
else
	R = sym(eye(n));
end
	c = cos(alpha_rad);
	s = sin(alpha_rad);
	R(idx,idx) =  c;
	R(jdx,jdx) =  c;
	R(idx,jdx) = -s; % signs were swapped
	R(jdx,idx) = +s;
end

