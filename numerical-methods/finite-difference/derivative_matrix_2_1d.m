% Sat Oct  2 21:21:48 MSD 2010
% Karl Kästner
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
% 
%
function D2 = derivative_matrix_2_1d(nx,L,order,bcl,bcr)
	if (nargin()<2||isempty(L))
		% choose domain [0 1]
		L = 1;
	end
	if (nargin() < 3 || isempty(order))
		order = 2;
	end
	if (nargin() < 4)
		bcl = 'dirichlet';
	end
	if (nargin() < 5)
		bcr = 'dirichlet';
	end
	if (isscalar(nx))
		% 
		n = nx;
	
		% note that this depends on the bnd:
		% for dirichlet only interior points
		% for neumann also exterior points
		%h = L/(n+1); % n? 1/(n-1); changed Sat Nov 27 02:01:15 MSK 2010
		h = L/(n-1);
		switch (order)
		case {4}
			% TODO there might be a mistake in the msc-thesis,
			%      as this kernel seems only valid for h=1
			D2 = [-1, 16, -30, 16, -1]/(12*h^2);
			m = 2;
		otherwise % == 2
			% second order: ddu_m := (1/h^2)*(u_l -2u_m + u_n);
			D2 = (1/h^2)*[1, -2, 1];
			m = 1;
		end
		D2 = spdiags(repmat(D2,n,1),-m:m,n,n);
		
		% boundaries
		switch (order)
		case {5}
			D2(1:2,:) = 0;
			D2(end-1:end,:) = 0;
			D2(2,1:3)       = 1/h^2*[1, -2, 1];
			D2(n-1,n-2:n)   = 1/h^2*[1, -2, 1];
			D2(1,1:3)       = 1/h^2*[1, -2, 1];
			D2(n,n-2:n)     = 1/h^2*[1, -2, 1];
		otherwise % 2
			D2(1,1:3)       = 1/h^2*[1, -2, 1];
			D2(n,n-2:n)     = 1/h^2*[1, -2, 1];
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
	end

if (0)
	switch (lower(bcl)) % homogeneous dirichlet f(0)=f(L)=0
		case {'dirichlet'}
			% nothing to do
		case {'neumann'} % homogeneous neumann df/dx(0) = 0, df/dx(L)=0
			D2(1,1:3) = 1/h^2*[1 -2 1];
			%2(1,1:4)       = 1/h^2*[ 2   -5    4   -1];
			%2(1,end-3:end) = 1/h^2*[-1    4   -5    2];
			%D2(1,1:5)       = 1/(12*h^2)*[35 -104  114  -56   11];
			%D2(1,end-4:end) = 1/(12*h^2)*[11  -56  114 -104   35];
		otherwise
			error('here');
	end % switch bcl
	switch (lower(bcr))
		case {'dirichlet'}
			% nothing to do
		case {'neumann'}
			D2(end,end-2:end) = 1/h^2*[1 -2 1];
		otherwise
			error('here');
	end% switch bcr
end
end % func derivative_matrix_2_1d

