% Fri 21 Aug 13:39:08 +08 2020
% central difference for scheme where values are kept at cell centres
%
	% central differences :
	% boundary condition
	% the boundary value is given at the cell interface,
	% but the variables at the cell centres
	% internal :       [  ..., -1/2,   0, 1/2, ... ]
	% series has to be teleskoping and thus conservative and second order at bc
	% the boundary difference kernel is thus (!)
	% dirichlet :
	% right :          [  ...,    0,-1/2, -1/2 ] + rval
	% left  :  -lval + [ +1/2, +1/2,   0,  ... ]
	% when the bc is f'|b_c = 0
	% extrap :
	%	lval = [3/2, -1/2, 0, ...]
	%	rval = [0, ..., -1/2, 3/2]
function dy = central_difference(x,y)
	% TODO heed variable step width
%	dy = [   -y(1) + 0.5*(y(2) + y(3));
%	        0.5*diff(y(3:end)-y(1:end-2));
%	      y(end) - 0.5*(y(end-1) + y(end-2))];
	dy = 0.5*(y(3:end) - y(1:end-2));
end
