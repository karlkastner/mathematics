% Sa 28. Nov 11:24:59 CET 2015
% Karl Kastner, Berlin
%% edges of a polygon
function E = poly_edges(xp,yp)
	np = length(xp);
	% edges of the polygon
	fdx = [0; find(isnan(xp)); np+1];
	E = [];
	for idx=1:length(fdx)-1
		n1 = fdx(idx)+1;
		n2 = fdx(idx+1)-1;
		% edges of the polygon and closing edge
		% TODO, this silently assumes the polygon was not yet closed
		if (n2 > n1)
			E = [E; [(n1:n2)', [(n1+1):n2,n1]']];
		end
	end
end

