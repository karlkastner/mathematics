%  Fr 13. Nov 17:28:57 CET 2015
%
%% centroid pf a polygone
function [cx, cy] = centroid(x,y)
	A  = Geometry.poly_area(x,y);
	cx = (1./(6*A)).*sum((x(1:end-1,:)+x(2:end,:)).*(x(1:end-1,:).*y(2:end,:) - x(2:end,:).*y(1:end-1,:)));
	cy = (1./(6*A)).*sum((y(1:end-1,:)+y(2:end,:)).*(x(1:end-1,:).*y(2:end,:) - x(2:end,:).*y(1:end-1,:)));
end

