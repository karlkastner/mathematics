% Apr 27 14:09
% Karl Kästner, Berlin 

function A = assemble_1d_dphi_dphi(P, T, func, int)
	nt1 = size(T,1);
	nt2 = size(T,2);
	nb = 0;

	% get quadrature coordinates and weigths
	[w b flag] = feval(int);

	% allocate memory
	% TODO, this is not enough for higher order
	buf = zeros(3*nt1,3);

	% calculate element integrals
	for idx=1:size(T,1)
		% fetch element coordinates
		A = [1 P(T(idx,1),:);
	             1 P(T(idx,2),:)];

		% determine the element length
		len = abs(det(A));

		% quadrature points in homogeneous Cartesian coordinates
		q = b*A;

		% construct the element point Vandermonde matrix spanning the test function polynomials
		% structure: 1 x x^2 x^3
		Va = vander_1d(P(T(idx,:)),nt2);

		% get test function coefficients
		C = inv(Va);

		% get coefficients of test function derivatives
		dC = derivative_1d(C, nt2-1);
		
		% evaluate test function derivative at quadrature points
		Vq = vander_1d(q(:,2), nt2-1);
		dphi = (Vq*dC);

		% evaluate the coefficient function at integration points
		if (nargin() > 2 && ~isempty(func))
			f = feval(func,q(:,2));
			% premultiply values
			wfa = w.*f*len;
		else
			% premultiply values
			wfa = w*len;
		end
		
		% add share to each of the elements corner points
		for adx=1:nt2
			% diagonal element
			I = -0.5*((wfa.*dphi(:,adx))'*dphi(:,adx));
			nb = nb+1;
			buf(nb,:) = [T(idx,adx) T(idx,adx) I];
			% off diagonal elements
			for bdx=adx+1:nt2
				% int (a1 + a2x)'*(b1 *b2x)'dx = int a2 b2 dx = a2 b2 * hx
				I = -((wfa.*dphi(:,adx))'*dphi(:,bdx));
				nb = nb+1;
				buf(nb,:) = [T(idx,adx) T(idx,bdx) I];
			end % for adx
		end % for bdx
	end % for idx
	% build matrix
	A = sparse(buf(:,1), buf(:,2),buf(:,3), size(P,1), size(P,1));
	% complete symmetric part
	A = A+A';
end % assemble_1d_dphi_dphi

