% Mo 9. Nov 11:26:22 CET 2015
% Karl Kastner, Berlin
%
%% lagrangian basis for interpolation
%% count number of valid samples
%
function [c, serr] = lp_count(x0,x,val)
	% TODO : so far order has to be constant or linear
	n   = length(x0);
	dof = n+1;

	% filter
	fdx = isfinite(val);
	x   = x(fdx);
	val = val(fdx);

%	% TODO this should be combined with the QR factorisation
%	% the current QR implementation fails for l(x) < dof
	if (length(x)<dof)
		c    = NaN(dof,1,class(val));
		serr = NaN(class(val));
		return;
	end

	% determine segment for each source point
	sdx = sum(bsxfun(@gt,cvec(x),rvec(x0)),2);
	sdx = min(n,max(1,sdx));

	c = zeros(n,0);
	for idx=1:n
		c(idx) = sum(sdx == idx);
	end
end % lp_regress

