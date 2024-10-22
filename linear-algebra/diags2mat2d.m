% 2024-10-03 10:46:04.561800380 +0200
% Karl Kastner, Berlin
%
%% construct sparse matrix from diagonals
%% similar to spdiags but does not shift rows
%% intended for circular boundary conditions
function A  = diags2mat2d(diags,n)
	switch (size(diags,2))
	case {5}
	% order u, l, r, d
	id  = reshape(1:prod(n),n);
	row = repmat(id(:),5,1);
	col = [flat(up(id,1)),flat(left(id,1)),flat(right(id,1)),flat(down(id,1)),flat(id)];
	A   = sparse(row(:),col(:),diags(:),n(1)*n(2),n(1)*n(2));
	case {9}
	% order clockwise
	% u, ur, r, dr, d, dl, l, ul
	id  = reshape(1:prod(n),n);
	row = repmat(id(:),9,1);
	col = [  flat(up(id,1)) ...
	       , flat(right(up(id,1),1)) ...
	       , flat(right(id,1)) ...
	       , flat(right(down(id,1),1)) ...
	       , flat(down(id,1)) ...
	       , flat(left(down(id,1),1)) ...
 	       , flat(left(id,1)) ...
	       , flat(left(up(id,1),1)) ...
	       , flat(id) ...
	      ];
	A   = sparse(row(:),col(:),diags(:),n(1)*n(2),n(1)*n(2));
	otherwise
		error('here');
	end
end

