% Tue May  1 01:19:46 MSK 2012
% Karl Kästner, Berlin
%
% integrates product of the derivatives of the test functions
%
% T : points on the triangle point, first three points are corners (size = [1 .. 6])
function [A buf] = assemble_2d_dphi_dphi(mesh, func, int)
	P = mesh.P;
	T = mesh.T;

	nt1 = size(T,1);
	nt2 = size(T,2);

	% degree of basis function polynomial (always an integer)
	% 1,3,6,10,15,21 -> 1,2,3,4,5,6
	%nv = -1.5 + sqrt(2*nt2 + 0.25);
	nv = (-3 + sqrt(8*nt2 + 1))/2;

	% get quadrature coordinates and weights
	% load the integration sample point coordinates for standard triangle
	[w b flag] = feval(int);

	% number of buffer entries per triangle
	mb0 = nt2*(nt2+1)/2;

	% exact number of elements in buffer after integration
	nb0 = nt1*mb0;

	% preallocate buffer
	buf = zeros(nb0,3);


	% integrate over each triangle
	nb = 0;
	for idx=1:nt1
		% corner points
		A = [   1 P(T(idx,1),:);
                        1 P(T(idx,2),:);
                        1 P(T(idx,3),:) ]; 

		% determine the triangle area
		area = 0.5*abs(det(A));

		% integration points in homogeneous Cartesian coordinates
		q = b*A;

		% point vandermonde matrix
		% structure: [1 x y xy x^2 y^2]
		Va = vander_2d(P(T(idx,:),:), nv);

		% compute test functions
		% C(i,:) * A(:,i) = 1; C(i,:) * A(:,j) = 0, i<>j
		% rows: test functions, columns: coefficients
		C = inv(Va);
	
		% form first derivative of the test functions
		% structure: dphi([1 x y]') : [dc00 dc01 dc10]*[1 x y]'
		[dC_x dC_y] = derivative_2d(C, nv);
	
		% evaluate test function derivative at the integration points
		% test function derivative evaluated at the points
		% rows: function, columns: points
		Vq = vander_2d(q(:,2:3), nv-1);
		dphi_x = Vq*dC_x;
		dphi_y = Vq*dC_y;

		% evaluate the coefficient function at the integration points
		% and premultiply values
		if (nargin() > 2 && ~isempty(func))
			wfa = -area*w.*feval(func,q(:,2:3));
		else
			wfa = -area*w;
		end

		% mass matrix contributions
		% for all testfunction being 1 at point adx
%		nb = mb0*(idx-1);
		for adx=1:nt2
			% diagonal entry integral approximation
			I = 0.5*sum( wfa.*(dphi_x(:,adx).^2 + dphi_y(:,adx).^2) );
			nb = nb+1;
			buf(nb,1) = T(idx,adx);
			buf(nb,2) = T(idx,adx);
			buf(nb,3) = I; % separate, as otherwise silentconversion to int 
			% off diagonal entries
			% exploit symmetry A(i,j) = A(j,i)
			for bdx=(adx+1):nt2
				% integral approximation
				I = sum( wfa.*(dphi_x(:,adx).*dphi_x(:,bdx) ...
						+ dphi_y(:,adx).*dphi_y(:,bdx)) );
				nb = nb+1;
				buf(nb,1) = T(idx,adx);
				buf(nb,2) = T(idx,bdx);
				buf(nb,3) = I; % separate, as otherwise silentconversion to int 
			end % for jdx
		end % for idx
	end % for idx (each triangle)

	if (nb ~= nb0)
		warning('preallocation insufficient');
	end

	% construct matrix from buffer
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(P,1), size(P,1));

	% complete symmetry
	A = A + A';
end % function assemble_2d_dphi_dphi

