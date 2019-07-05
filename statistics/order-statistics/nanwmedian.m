% Di 7. Apr 16:41:58 CEST 2015
% Karl Kastner, Berlin
%
%% weighted median, skips nan-values
function [me s l u] = nanwmedian(w,x,P)
	if (nargin() < 3 || isempty(P))
		% P = 0.95;
		% one-sigma interval
		P = normcdf(1);
	end

	if (isvector(w))
		w = repmat(cvec(w),1,size(x,2));
	end

	%w(~isfinite(x)) = NaN;
	%n = nansum(w);
	if (size(w,1) ~= size(x,1))
		error('w and x must be of same length');
	end

	n = max(1,size(x,2));
	me = NaN(1,n);
	if (nargout()<2)
		for idx=1:size(x,2)
			me(:,idx) = nanwquantile(w(:,idx),x(:,idx),0.5);
		end
	else
		dof = wdof(w);
		[pl pm pu] = median_ci(dof,P);
		for idx=1:size(x,2)
			q(:,idx) = nanwquantile(w(:,idx),x(:,idx),[pl pm pu]);
		end
		me   = q(2,:);
		l    = q(1,:);
		u    = q(3,:);
		s    = 0.5*(u-l);
	end
end

