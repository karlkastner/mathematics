% 2015-06-02 09:18:59.390325255 +0200
% Karl Kastner, Berlin
% analytic histogram
function obj = setup(obj,cdf)
	edges    = obj.edge;
	H	 = cdf(edges);
	n        = size(H,1);
	% extend to infinity to make integral 1
	H(:,1)   = 0;
	H(:,end) = 1;
%	H	 = [zeros(n,1), H, ones(n,1)];
	obj.h        = diff(H,[],2);
%	h(:,1)   = cdf(edges(2));
%	for idx=2:length(edges)-1
%		h(:,idx) = cdf(edges(idx+1))-cdf(edges(idx));
%	end
%	h(:,end) = 1-cdf(edges(end-1));
end

