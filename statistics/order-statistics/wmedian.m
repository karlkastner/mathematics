% Di 7. Apr 16:41:58 CEST 2015
% Karl Kastner, Berlin
%
%% weighted median
function [me, s, l, u] = wmedian(w,x,P)
	if (nargin() < 3 || isempty(P))
		% P = 0.95;
		% one-sigma interval
		P = normcdf(1);
	end

	if (isvector(w))
		w = repmat(cvec(w),1,size(x,2));
	end
	if (size(w,1) ~= size(x,1))
		error('w and x must be of same length');
	end
	
	me = NaN(1,size(x,2));
	if (nargout()<2)
		for idx=1:size(x,2)
			me(:,idx) = wquantile(w(:,idx),x(:,idx),0.5);
		end
	else
		dof = wdof(w);
		[pl, pm, pu] = median_ci(dof,P);
		q   = zeros(3,size(x,2));
		for idx=1:size(x,2)
			q(:,idx) = wquantile(w(:,idx),x(:,idx),[pl, pm, pu]);
		end
		me   = q(2,:);
		l    = q(1,:);
		u    = q(3,:);
		% TODO this is only true for P = normcdf(1)
		s    = 0.5*(u-l);
	end
end

