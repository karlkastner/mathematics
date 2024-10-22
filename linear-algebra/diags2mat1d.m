% Fri  4 Oct 11:41:21 CEST 2024
% Karl Kastner, Berlin
%% set up tridiagonal matrix, intended for circular boundary conditions
%% similar to spdiags but rows are not shifted
function A  = diags2mat1d(diags,n)
	id  = (1:n)';
	row = repmat(id(:),3,1);
	col = [flat(up(id,1)),flat(down(id,1)),flat(id)];
	A   = sparse(row(:),col(:),diags(:),n,n);
end

