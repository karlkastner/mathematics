% Tue 25 Oct 17:19:07 CEST 2016
% Karl Kastner, Berlin
%% height of a triangle
function h = tri_height(x,y);
	l = Geometry.tri_edge_length(x,y);
	s = sum(l,2)/2;
	h_ = 2*sqrt(s.*(s-l(:,1)).*(s-l(:,2)).*(s-l(:,3)));
	h = bsxfun(@times,h_,1./l);
end

