% Mon 23 Jul 11:18:08 CEST 2018
function [p, obj] = cdf(obj,q)
	H    = obj.cdfS(obj.h);
	p    = zeros(size(obj.h,1),length(q));
	cdf = obj.cdf();
	for idx=1:size(obj.h,1)
		p(idx,:) = interp1(cdf(:,idx),obj.centre,p,'linear','extrap');
	end
	p(p<0) = 0;
	p(p>1) = 1;
end

