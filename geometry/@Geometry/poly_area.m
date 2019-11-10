% Fr 13. Nov 17:28:45 CET 2015
%% area of a polygon
%% function A = poly_area(x,y)
function A = poly_area(x,y)
	% shift for better conditioning
	x = bsxfun(@minus,x,x(1,:));
	y = bsxfun(@minus,y,y(1,:));

	A = 0.5*sum(x(1:end-1,:).*y(2:end,:) - x(2:end,:).*y(1:end-1,:));
end

