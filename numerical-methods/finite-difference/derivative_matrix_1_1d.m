% Sat Oct  2 21:12:57 MSD 2010
% changed  Sat Nov 27 04:39:21 MSK 2010
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%% finite difference matrix of first derivative in one dimensions
%% n : number of grid points
%% h = L/(n+1) constant step with
%%
%% function [D1, d1] = derivative_matrix_1_1d(arg1,arg2,order,bc,bcr,isdx)
function [D1, d1] = derivative_matrix_1_1d(arg1,arg2,order,bcl,bcr,isdx)
	if (nargin()<6)
		isdx = false;
	end
	if (nargin()< 2 || isempty(arg2))
		% choose domain [0 1]
		arg2 = 1;
	end

	if (nargin() < 3 || isempty(order))
		order = 2;
	end
	if (nargin() < 4 || isempty(bcl))
		bcl = 'hdirichlet';
	end
	if (nargin() < 5 || isempty(bcr))
		bcr = bcl;
	end
	switch(bcl)
	case {'neumann','hdirichlet','circular'}
		% nothing to do
	otherwise
		error('unknown boundary condition');
	end

	if (isscalar(arg1))

	% equal grid spacing
	n = arg1;

	% changed Sat Nov 27 02:01:15 MSK 2010
	if (~isdx)
		L = arg2;
		h = L / (n-1);
	else
		h = arg2;
	end


	switch (order)
	 case {-1,'-1'}
		% first order backward
		% du[i] = (1/h)(u[i+1] - u[i])
		d1 = 1.0/h*[0, -1,  1];
	 case {+1,'+1'}
		% first order forward
		% du[i] = (1/h)(u[i-1] - u[i])
		d1 = 1.0/h*[-1,  1,  0];
	case {'+2'}
		% second order left
		d1 = 1.0/h*[1/2, -2, 3/2];
	case {2,'2'}
		% second order central
		% du[i] = (1/2h)*(u[i+1] - u[i-1])
		d1 = 0.5/h*[-1,  0,  1];
	case {-2,'-2'}
		% second order right
		d1 = 1.0/h*[-3/2, 2, -1/2];
	case {3,'3',4}
		d1 = 1/(12*h)*[1, -8, 0, 8, -1];
	case {'+3'}
		k = -2:0;
		d1=1/h*[ -1,   18,   63,   46]/30;
	case {'-3'}
		k = 0:3;
		d1=1/h*[-46,   63,  -18,    1]/30;
	case {6}
		d1 = 1/(60*h)*[-1, 9, -45, 0, 45, -9, 1];
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
%		if (~strcmp(bc,'circular'))
%			d1(end,2)   =   1/h;
%			d1(end-1,1) =  -1/h;
%		end
		D1 = spdiags(d1,-1:1,n,n);
		switch (bcl)
		case ('hdirichlet')
			% nothing to do, outside value is zero
			% TODO this depends where the boundary is placed
		case ('neumann')
			% [0,-1| 1]/h	bc and pde
			% [-1,2|-1]/h	extrapolation
			% [-1,1]/h	
			D1(end,end-1:end) = [-1,1]/h;
		case ('circular')
			D1(end,1) = 1/h;
		end
	case {+1,'+1'}
%		if (~strcmp(bc,'circular'))
%			d1(1,2) = -1/h;
%			d1(2,3) =  1/h;
%		end
		D1 = spdiags(d1,-1:1,n,n);
%		if (strcmp(bc,'circular'))
%			D1(1,end) = -1/h;
%		end
		switch (bcr)
		case ('hdirichlet')
			% nothing to do, outside value is zero
			% TODO this depends where the boundary is placed
		case ('neumann')
			% [-1| 1, 0]/h	bc and pde
			% [-1| 2,-1]/h	extrapolation
			% [0 |-1, 1]/h
			D1(1,1:2) = [-1,1]/h;
		case ('circular')
			D1(1,end) = -1/h;
		end
	case {'+2'}
		D1 = spdiags(d1,-2:0,n,n);
		D1(1,1:2) = 1/h*[-1, 1];
		D1(2,1:2) = 1/h*[-1, 1];
	case {2,'2'}
		if (~strcmp(bcl,'circular'))
			d1(1,2)     = -1/h;
			d1(2,3)     =  1/h;
		end	
		if (~strcmp(bcr,'circular'))
			d1(end-1,1) = -1/h;
			d1(end,2)   =  1/h;
		end
		D1 = spdiags(d1,-1:1,n,n);
		if (strcmp(bcl,'circular'))
			D1(1,end) = -0.5/h;
			D1(end,1) = +0.5/h;
		end
		% note there seems to be a mistake in table of my master thesis!
		%D1(1,1:3) = 0.5/h*[-3 4 -1];
		%D1(end,end-2:end) = 0.5/h*[1 -4 3];
	case {-2,'-2'}
		D1 = spdiags(d1,0:2,n,n);
		D1(end-1,end-1:end) = 1/h*[-1, 1];
		D1(end,end-1:end)   = 1/h*[-1, 1];
	case {'+3','-3'}
		D1 = spdiags(d,k,n,n);
	case {3,'3',4}
		D1 = spdiags(d1,-2:2,n,n);

		if (strcmp(bcr,'circular'))
				D1(2,end)   = -D1(2,4);
				D1(end-1,1) = -D1(end-1,end-3);

				D1(1,end)   = -D1(1,2);
				D1(1,end-1) = -D1(1,3);
				D1(end,2)   = -D1(end,end-2);		
				D1(end,1)   = -D1(end,end-1);	
		else
		warning('bc not yet implemented')
		end
	case {6}
		D1 = spdiags(d1,-3:3,n,n);
		if (strcmp(bcl,'circular'))
				D1(3,end)   = -D1(3,6);
				D1(end-2,1) = -D1(end-2,end-5);

				D1(2,end)    = -D1(2,4);
				D1(2,end-1)  = -D1(2,5);
				D1(end-1,2)  = -D1(end-1,end-4);		
				D1(end-1,1)  = -D1(end-1,end-3);
			
				D1(1,end)    = -D1(1,2);
				D1(1,end-1)  = -D1(1,3);
				D1(1,end-2)  = -D1(1,4);
				D1(end,3)    = -D1(end,end-3);		
				D1(end,2)    = -D1(end,end-2);		
				D1(end,1)    = -D1(end,end-1);
		else
			warning('bc not yet implemented')
		end
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

