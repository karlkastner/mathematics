% 15:05:40.873225890 +0200
% Karl Kastner, Berlin
%
%% extract 5-point diagonals from 2D discretization matrix with poten
function diags = extract_diagonals(A,n)
	id = reshape(1:n(1)*n(2),n(1),n(2));
	
	u = up(id,1);
	d = down(id,1);
	l = left(id,1);
	r = right(id,1);

	id = sub2ind(n(1)*n(2)*[1,1],repmat(id(:),5,1), [u(:);l(:);id(:);r(:);d(:)]);

	diags = reshape(full(A(id)),n(1)*n(2),5);
end

