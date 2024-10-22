% Tue 17 Sep 14:58:33 CEST 2024
% Karl Kastner, Berlin
%
%% assemble matrix from colums of a 1D 3-point kernel
%% similar to spdiags, but does not shift colums
%% predecessor of first and successor of last elements are wrapped
%% (for circular boundary conditions)
function A = diags2mat_1d(diags,n)
	id = (1:n)';
	
	u = up(id,1);
	d = down(id,1);

	row = repmat(id(:),3,1);
	col = [u(:);id(:);d(:)];
	A = sparse(row,col,diags(:),n,n);
end

