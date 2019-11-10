% Fri  8 Nov 16:36:43 +08 2019
function cflag = isconvex_quad(X,Y)
	A = zeros(size(X));
	% the four sub-triangles have the same area sign
	for idx=1:4
		A(:,idx) = Geometry.tri_area(X(:,1:3),Y(:,1:3));
		X = [X(:,2:4),X(:,1)];
		Y = [Y(:,2:4),Y(:,1)];
	end
	p     = prod(A,2);
	cflag = p>0;
end

