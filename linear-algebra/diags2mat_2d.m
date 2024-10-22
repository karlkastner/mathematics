% Tue 17 Sep 14:58:33 CEST 2024
% Karl Kastner, Berlin
%
%% assemble matrix from colums of a 2D 5-point kernel
%% similar to spdiags, but does not shift colums
%% predecessor of first and successor of last elements are wrapped
%% (for circular boundary conditions)
function A = diags2mat_2d(diags,n)
	nn = n(1)*n(2);
	id = reshape(1:nn,n(1),n(2));
	
	u = up(id,1);
	d = down(id,1);
	l = left(id,1);
	r = right(id,1);

	%id = sub2ind(n(1)*n(2)*[1,1],repmat(id(:),5,1), [u(:);l(:);id(:);r(:);d(:)]);
	row = repmat(id(:),5,1);
	col = [u(:);l(:);id(:);r(:);d(:)];
	A = sparse(row,col,diags(:),nn,nn);
end	

