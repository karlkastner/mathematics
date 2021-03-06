% Sat Oct  2 21:12:57 MSD 2010
% changed  Sat Nov 27 04:39:21 MSK 2010
% Karl Kästner, Berlin
%
%% finite difference matrix of first derivative in one dimensions
%% n : number of grid points
%% h = L/(n+1) constant step with
%% function [D1, d1] = derivative_matrix_1d(n,L,order)
% TODO allow optionally for circular boundary condition
function [D1, d1] = derivative_matrix_1_1d(arg1,L,order)
	if (nargin()< 2 || isempty(L))
		% choose domain [0 1]
		L = 1;
	end
	if (nargin() < 3 || isempty(order))
		order = 2;
	end

	if (isscalar(arg1))

	% equal grid spacing
	n = arg1;

	% changed Sat Nov 27 02:01:15 MSK 2010
	%h = L/(n+1); % n? 1/(n-1);
	h = L/(n-1);

	switch order
	 case {-1,'-1'}
		% first order backward
		% du[i] = (1/h)(u[i+1] - u[i])
		d1 = 1.0/h*[0 -1  1];
	 case {+1,'+1'}
		% first order forward
		% du[i] = (1/h)(u[i-1] - u[i])
		d1 = 1.0/h*[-1  1  0];
	case {'+2'}
		% second order left
		d1 = 1.0/h*[1/2, -2, 3/2];
	case {2,'2'}
		% second order central
		% du[i] = (1/2h)*(u[i+1] - u[i-1])
		d1 = 0.5/h*[-1  0  1];
	case {-2,'-2'}
		% second order right
		d1 = 1.0/h*[-3/2, 2, -1/2];
	 case {3,'3'}
		d1 = 1/(12*h)*[1 -8 0 8 -1];
	 otherwise
		error('not yet implemented');
	end % switch

	% repeat difference kernel for each row
	d1 = repmat(d1,n,1);

	% diagonals to matrix
	% construct matrix from diagonals
	% and correct first and last rows at boundary
	switch (order)
	case {-1,'-1'}
		d1(end,2)   =   1/h;
		d1(end-1,1) =  -1/h;
		D1 = spdiags(d1,-1:1,n,n);
	case {+1,'+1'}
		d1(1,2) = -1/h;
		d1(2,3) =  1/h;
		D1 = spdiags(d1,-1:1,n,n);
	case {'+2'}
		D1 = spdiags(d1,-2:0,n,n);
		D1(1,1:2) = 1/h*[-1, 1];
		D1(2,1:2) = 1/h*[-1, 1];
	case {2,'2'}
		d1(1,2)     = -1/h;
		d1(2,3)     =  1/h;
		d1(end-1,1) = -1/h;
		d1(end,2)   =  1/h;
		D1 = spdiags(d1,-1:1,n,n);
		% note there seems to be a mistake in table of my master thesis!
		%D1(1,1:3) = 0.5/h*[-3 4 -1];
		%D1(end,end-2:end) = 0.5/h*[1 -4 3];
	case {-2,'-2'}
		D1 = spdiags(d1,0:2,n,n);
		D1(end-1,end-1:end) = 1/h*[-1, 1];
		D1(end,end-1:end)   = 1/h*[-1, 1];
	case {3,'3'}
		D1 = spdiags(d1,-2:2,n,n);
	otherwise
		error('here');
	end

	else	
		x = cvec(arg1);
		n = length(x);
		% grid spacing
		h = diff(x);

		switch (order)
		case {-1,'-1'}
			d1 = [[-1./h;0],[0;1./h]];
			D1 = spdiags(d1,-1:0,n,n);
			D1(1,1:2) = [-1./h(1),1./h(1)];
		case {+1,'+1'}
			%d1 = [-h,h;-h(end),h(end)];
			d1 = [[-1./h;0],[0;1./h]];
			D1 = spdiags(d1,0:1,n,n);
			D1(end,end-1:end) = [-1./h(end),1./h(end)];
		case {2,'2'}
			hr = [0;h(2:end);0];
			hl = [0;h(1:end-1);0];
			d1 = [ -hr./(hl.*(hl+hr)), ...
			       -(hl-hr)./(hl.*hr), ...
			       +(hl)./(hr.*(hl+hr)); ...
			];
 	               D1 = spdiags( [ [d1(2:end,1); 0] d1(:,2) [0; d1(1:end-1,3)]], -1:1, n, n);
		       % TODO, these are just first order accurate
		       D1(1,1:2) = [-1./h(1),1./h(1)];
		       D1(end,end-1:end) = [-1./h(end),1./h(end)];
		otherwise
			error('not yet implemented');
		end
	end
end % function derivative_matrix_1_1d

