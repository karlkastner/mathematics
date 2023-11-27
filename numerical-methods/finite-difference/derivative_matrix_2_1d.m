% Sat Oct  2 21:21:48 MSD 2010
% Karl Kästner
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% finite derivative matrix of second derivative in one dimension
% 
%
% boundaries need to be taken care of by user
% only second order supported for second derivate
%
% TODO use derive_kernel to set up matrices on the fly
% TODO it would be best not to eliminate the first and last point for dirichlet bc, to have consistent matrix sizes
%
% nx : scalar : number of grid points or grip points sorted, for equispaced grid
% nx : vector : sorted grid points, not necessarily sorted
% L  : lenght of domain, will not be used if nx is scalar
% order : [2,4] : accuracy of discretisation, only second order for variable grid 
% bcl : left boundary condition
%	hdirichlet : homogeneous dirichlet   y(0) = 0
%	idirichlet : inhomogeneous dirichlet y(0) = y0
%	neumann    : homogeneous neumann     dy/dx(0) = 0
%	TODO	     inhomogeneous neumann   dy/dx(0) = dydx0
%	circular   : circular                y(0) = y(L+h)
% bcr : right boundary condition, defined as above 
%
% function D2 = derivative_matrix_2_1d(nx,L,order,bcl,bcr)
function D2 = derivative_matrix_2_1d(nx,arg2,order,bcl,bcr,isdx)
	if (nargin()<6)
		isdx = false;
	end
	if (nargin()<2||isempty(arg2))
		% choose domain [0 1]
		arg2 = 1;
	end
	if (~isdx)
		h = arg2/(nx-1);
	else
		h = arg2;
	end
	if (nargin() < 3 || isempty(order))
		order = 2;
	end
	if (nargin() < 4)
		% homogeneous dirichlet
		bcl = 'hdirichlet';
	end
	if (nargin() < 5)
		bcr = bcl;
	end
	if (isscalar(nx))
	
		% note that this depends on the bnd:
		% for dirichlet only interior points
		% for neumann also exterior points
		%h = L/(n+1); % n? 1/(n-1); changed Sat Nov 27 02:01:15 MSK 2010
		switch (order)
		case {2}
			% second order: ddu_m := (1/h^2)*(u_l -2u_m + u_n);
			D2 = (1/h^2)*[1, -2, 1];
			m = 1;
		case {4}
			% TODO there might be a mistake in the msc-thesis,
			%      as this kernel seems only valid for h=1
			D2 = [-1, 16, -30, 16, -1]/(12*h^2);
			m = 2;
		case {6}
			D2 = [2,   -27,   270,  -490,    270,   -27,     2]/(180*h*h);
			m = 3;	
		otherwise 
			error('order must be 2, 4 or 6');
		end
		D2 = spdiags(repmat(D2,nx,1),-m:m,nx,nx);

		% left boundary
		switch (lower(bcl))
		case {'hdirichlet'}
			switch (order)
			case {2}
				% nothing to do
			otherwise
				error('not yet implemented');
			end	
		case {'idirichlet'}
			switch (order)
			case {2}
				D2(1,:) = 0;
				D2(1,1) = 1;
			otherwise
				error('not yet implemented');
			end
		case {'circular'}
			switch (order)
			case {2}
				D2(1,end) = D2(1,2);
				%D2(end,1) = D2(end,end-1);
			case {4}
				D2(2,end)   = D2(2,4);
				%D2(end-1,1) = D2(end-1,end-3);

				D2(1,end)   = D2(1,2);
				D2(1,end-1) = D2(1,3);
				%D2(end,2)   = D2(end,end-2);		
				%D2(end,1)   = D2(end,end-1);	
			case {6}
				D2(3,end)   = D2(3,6);
				%D2(end-2,1) = D2(end-2,end-5);

				D2(2,end)    = D2(2,4);
				D2(2,end-1)  = D2(2,5);
				%D2(end-1,2)  = D2(end-1,end-4);		
				%D2(end-1,1)  = D2(end-1,end-3);
			
				D2(1,end)    = D2(1,2);
				D2(1,end-1)  = D2(1,3);
				D2(1,end-2)  = D2(1,4);
				%D2(end,3)    = D2(end,end-3);		
				%D2(end,2)    = D2(end,end-2);		
				%D2(end,1)    = D2(end,end-1);
			otherwise
				error('not yet implemented');
			end % switch order
		case {'neumann'}
			% linear extrapolation
			switch (order)
			case {2}
				%     [ 1 | -2, 1]/h^2
				% bc: [-1 | 1]/h == 0
				%     [ 0 | -1, 1]/h^2 = 0
				%D2(1,1:3)       = 1/h^2*[1, -2, 1];
				%D2(n,n-2:n)     = 1/h^2*[1, -2, 1];
				D2(1,1:2)        = 1/(h*h)*[-1,1];
				%    [1, -2 | 1]/h^2 (ode)
				%    [   -1 | 1] = 0 (bc)
				%    [1, -1 | 0]/h^2
				%D2(end,end-1:end) = 1/(h*h)*[1,-1]; % was -1, 1
			case {4}
				D2(1:2,:) = 0;
				%D2(end-1:end,:) = 0;
				D2(2,1:3)       = 1/h^2*[1, -2, 1];
				%D2(n-1,n-2:n)   = 1/h^2*[1, -2, 1];
				D2(1,1:3)       = 1/h^2*[1, -2, 1];
				%D2(n,n-2:n)     = 1/h^2*[1, -2, 1];
			otherwise
				error('not yet implemented');
			end % switch order
		otherwise
			error('not yet implemented');
		end

		% right boundary
		switch (lower(bcr))
		case {'hdirichlet'}
			switch (order)
			case {2}
				% nothing to do
			otherwise
				error('not yet implemented');
			end	
		case {'idirichlet'}
			switch (order)
			case {2}
				D2(nx,:) = 0;
				D2(nx,nx) = 1;
			otherwise
				error('not yet implemented');
			end
		case {'circular'}
			switch (order)
			case {2}
				D2(end,1) = D2(end,end-1);
			case {4}
				D2(end-1,1) = D2(end-1,end-3);

				D2(end,2)   = D2(end,end-2);		
				D2(end,1)   = D2(end,end-1);	
			case {6}
				D2(end-2,1) = D2(end-2,end-5);

				D2(end-1,2)  = D2(end-1,end-4);		
				D2(end-1,1)  = D2(end-1,end-3);
			
				D2(end,3)    = D2(end,end-3);		
				D2(end,2)    = D2(end,end-2);		
				D2(end,1)    = D2(end,end-1);
			end % switch order
		case {'neumann'}
			% linear extrapolation
			switch (order)
			case {2}
				%     [ 1 | -2, 1]/h^2
				% bc: [-1 | 1]/h == 0
				%     [ 0 | -1, 1]/h^2 = 0
				%D2(1,1:3)       = 1/h^2*[1, -2, 1];
				%D2(n,n-2:n)     = 1/h^2*[1, -2, 1];
				% D2(1,1:2)        = 1/(h*h)*[-1,1];
				%    [1, -2 | 1]/h^2 (ode)
				%    [   -1 | 1] = 0 (bc)
				%    [1, -1 | 0]/h^2
				D2(end,end-1:end) = 1/(h*h)*[1,-1]; % was -1, 1
			case {4}
				D2(nx-1:nx,:) = 0;
				D2(nx-1,nx-2:nx)   = 1/h^2*[1, -2, 1];
				D2(nx,nx-2:nx)     = 1/h^2*[1, -2, 1];
			otherwise
				error('not yet implemented');	
			end
		otherwise
			error('not yet implemented');
		end
	else
		x = cvec(nx);
		n = length(x);

            	xl = [0; x(1:end-1)];
	        xc = x;
        	xr = [x(2:end);0];

                A = 2./[ (xr-xl).*(xc-xl), ...                                  
                    -(xr-xc).*(xc-xl), ...                                      
                     (xr-xc).*(xr-xl) ];                                        
                D2 = spdiags( [ [A(2:end,1); 0] A(:,2) [0; A(1:end-1,3)]], -1:1, n, n);
		% TODO is this exact for variable meshes? check with thesis
		D2(1,1:3) = D2(2,1:3);
		D2(end,end-2:end) = D2(end-1,end-2:end);
			% TODO implement other boundary conditions
	end
end % func derivative_matrix_2_1d

