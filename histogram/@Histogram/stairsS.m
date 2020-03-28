function [wh_ e] = stairS(h,edge_A)
%	h = obj.h;
%	edge_A = obj.edge;
	n = size(h,1);
	m = size(h,2);
	% width normalisation
	w = abs(rvec(diff(edge_A)));
%	w = 1./mean(w).*w;
	wh = bsxfun(@times, 1./w, h);
	e = zeros(1,2*m);
	wh_ = zeros(n,2*m);
	for idx=1:m
		wh_(:,2*idx-1:2*idx,1) = wh(:,idx);
		e(2*idx-1) = edge_A(idx);
		e(2*idx)   = edge_A(idx+1);
	end
end

