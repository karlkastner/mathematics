% Sat Feb 25 23:51:16 MSK 2012
% Karl Kästner, Berlin

function A = assemble_3d_dphi_dphi(P, T, func, int)
	nt1 = size(T,1);
	nt2 = size(T,2);

	% degree of basis function polynomial (always an integer)
	% 1,4,10,20,35,... -> 0,1,2,3,4,...
	% lt2 = 1/6*nv^3 + nv^2 + 11/6*nv + 1
 	nv = round(1/(3*(3*nt2+(9*nt2^2-1/27)^(1/2))^(1/3))+(3*nt2+(9*nt2^2-1/27)^(1/2))^(1/3)-2);

	% get quadrature coordinates and weights
	[w b flag] = feval(int);

	% number of buffer elements per triangle
	mb0 = nt2*(nt2+1)/2;

	% exact number of elements in buffer after integration
	nb0 = nt1*mb0;

	% allocate memory
	buf = zeros(nb0,3);

	% calculate element integrals
	nb = 0;
	for idx=1:nt1
		% fetch element coordinates
		A = [1 P(T(idx,1),:);
	             1 P(T(idx,2),:);
	             1 P(T(idx,3),:);
		     1 P(T(idx,4),:)];

		% strange - laplacian eigs require 0.5, harm-osci requires 0.5^2
		volume = 1.0/6.0*abs(det(A));

		% quadrature points in homogeneous Cartesian coordinates
		q = b*A;

		% Vandermonde matrix of support points of the polynomial basis function
		Va = vander_3d(P(T(idx,:),:), nv);

		% get test function coefficients
		C = inv(Va);

		% get the derivative
		[dC_x dC_y dC_z] = derivative_3d(C);

		% evaluate the derivative at the quadrature points
		Vq = vander_3d(q(:,2:4),nv-1);
		dphi_x = Vq*dC_x;
		dphi_y = Vq*dC_y;
		dphi_z = Vq*dC_z;

		% evaluate coefficient function at quadrature points
		% and precompute weigths
		if (nargin() > 3 && ~isempty(func))
			wfa = -volume*w.*feval(func, q(:,2:4));
		else
			wfa = -volume*w;
		end

		% add share to each of the elements corner points
%		nb = mb0*(idx-1);
		for adx=1:nt2
			% diagonal entry
			I = 0.5*(         (wfa.*dphi_x(:,adx))'*dphi_x(:,adx) ...
					+ (wfa.*dphi_y(:,adx))'*dphi_y(:,adx) ...
					+ (wfa.*dphi_z(:,adx))'*dphi_z(:,adx) );
			nb = nb+1;
			%A(p(adx),p(bdx)) = A(p(adx),p(bdx)) + I;
			buf(nb,1) = T(idx,adx);
			buf(nb,2) = T(idx,adx);
			buf(nb,3) = I;
			% off diagonal entries
			for bdx=adx+1:nt2
				I = (     (wfa.*dphi_x(:,adx))'*dphi_x(:,bdx) ...
					+ (wfa.*dphi_y(:,adx))'*dphi_y(:,bdx) ...
					+ (wfa.*dphi_z(:,adx))'*dphi_z(:,bdx) );
				nb = nb+1;
				buf(nb,1) = T(idx,adx);
				buf(nb,2) = T(idx,bdx);
				buf(nb,3) = I;
			end % for adx
		end % for bdx
	end % idx
	
	if (nb ~= nb0)
		warning('assemple_3d_dphi_dphi', 'preallocation insufficient');
	end

	% build matrix
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(P,1), size(P,1));

	% symetric completition
	A = A + A';
end % assemble_3d_dphi_dphi

