% Tue May  1 01:20:41 MSK 2012
% Karl Kästner, Berlin
%
% TODO exact integration with f=1
%
% T : global point indices (size = [1 6])
% P : point coordinates (size)
% f : potential function (function handle)
function [A buf] = assemble_2d_phi_phi(mesh, func, int)
	P = mesh.P;
	T = mesh.T;

	nt1 = size(T,1);
	nt2 = size(T,2);

	% degree of basis function polynomial (always an integer)
%	nv = -1.5 + sqrt(2*nt2 + 0.25);
	nv = (-3 + sqrt(8*nt2 + 1))/2;

	% load the integration sample point coordinates for standard triangle
	[w b flag] = feval(int);

	if (0 ~= flag && isempty(func))
		% diagonal lumped mass matrix
		% test functions phi_i are chosen that integral is 1 at the points p_i and zero at the other points
		% number of buffer entries per triangle
		mb0 = nt2;
	else
		% non diagonal mass matrix
		% number of buffer entries per triangle
		mb0 = nt2*(nt2+1)/2;
	end

	% exact number of elements in buffer after integration
	nb0 = nt1*mb0;

	% preallocate buffer
	buf = zeros(nb0,3);

	% integrate over each triangle
	% embarassingly parallel and still matlab is too stupid to allow parfor here
	nb = 0;
	for idx=1:nt1
		% corner points
		A = [   1 P(T(idx,1),:);
                        1 P(T(idx,2),:);
                        1 P(T(idx,3),:) ]; 

		% determine the triangle area
		area = 0.5*abs(det(A(1:3,1:3)));

		% integration points in homogeneous Cartesian coordinates
		q = b*A;

		% construct the triangle point Vandermonde matrix spanning the test function polynomials
		% structure: 1 x y xy x^2 y^2 x^2y xy^2 x^3 y^3
		Va = vander_2d(P(T(idx,:),:), nv);

		% compute test functions
		% C(i,:) * A(:,i) = 1; C(i,:) * A(:,j) = 0, i<>j
		C = inv(Va);

		% evaluate test function at the integration points
		% TODO these values should be constant and can be precomputed
		Vq = vander_2d(q(:,2:3), nv);

		% phi is not perfectly symmetric!!!
		phi = Vq*C;

		% evaluate the coefficient function at integration points
		if (nargin() > 2 && ~isempty(func))
			f = feval(func,q(:,2:3));
			% premultiply values
			wfa = w.*f*area;
		else
			% premultiply values
			wfa = w.*area;
		end
		
		% stiffness matrix contribution
		% for all testfunction being 1 at point adx
		nb = mb0*(idx-1);
		for adx=1:nt2
			% diagonal entry integral approximation
			I = 0.5*sum(wfa.*phi(:,adx).^2);
			nb = nb+1;
			buf(nb,1) = T(idx,adx);
			buf(nb,2) = T(idx,adx);
			buf(nb,3) = I; % not in common array, otherwise silent conversion into int
			% off diagonal entries
			% exploit symmetry A(i,j) = A(j,i)
			if (0 == flag || 0 == isempty(func))
				for bdx=(adx+1):nt2
					% integral approximation
					I = sum(wfa.*phi(:,adx).*phi(:,bdx));
					nb = nb+1;
					buf(nb,1) = T(idx,adx);
					buf(nb,2) = T(idx,bdx);
					buf(nb,3) = I; % not in common array, otherwise silent conversion into int
				end % for jdx
			end
		end % for idx
	end % for idx (each triangle)

	if (nb ~= nb0)
		warning('preallocation insufficient');
	end
	% construct the matrix from the buffer
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(P,1), size(P,1));
	% complete symmetry
	A = A + A';
end % assemble_2d_phi_phi

