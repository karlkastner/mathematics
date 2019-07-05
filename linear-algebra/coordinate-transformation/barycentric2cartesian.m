% Wed  8 Jun 09:25:16 CEST 2016
%
%% barycentric to cartesian coordinates
function [xy] = barycentric2cartesian(xy0,pq)
	pq = pq(:,1:2);
	pq = pq';
%	pq   = reshape(pq,2,[])
	T    = [pq; 1-pq(1,:)-pq(2,:)]';
%	xy0_ = reshape(xy0,2,[])'
	xy0_ = xy0;
	xy   = T*xy0_;
	xy   = xy';
	%xy   = xy(:);
end


