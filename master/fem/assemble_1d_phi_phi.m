% Apr 27 13:54 MSD
% Karl Kästner, Berlin

function A = assemble_1d_phi_phi(P, T, func, int)
	nt = size(T,1);
	nt2 = size(T,2);
	nb = 0;

	% get quadrature coordinates and weigths
	[w b flag] = feval(int);
	
	% allocate memory
	% TODO, this is not enough for higher order
	if (flag)	
		buf = zeros(2*nt,3);
	else
		buf = zeros(3*nt,3);
	end

	% integrate shares of each triangle
	for idx=1:size(T,1)
		% fetch element coordinates
		A = [   1 P(T(idx,1));
                        1 P(T(idx,2)); ];

		% determine the element length
		len = abs(det(A));

		% integration points in homogeneous Cartesian coordinates
		q = b*A;

		% construct the element point Vandermonde matrix spanning the test function polynomials
		% structure: 1 x x^2 x^3
		Va = vander_1d(P(T(idx,:)),nt2);

		% get test function coefficients
		C = inv(Va);

		% evaluate test function at the integration points
		% TODO these values should be constant and can be precomputed
		Vq = vander_1d(q(:,2), nt2);

		% phi is not perfectly symmetric!!!
		phi = (Vq*C);

		% evaluate the coefficient function at integration points
		if (nargin() > 2 && ~isempty(func))
			f = feval(func,q(:,2));
			% premultiply values
			wfa = w.*f*len;
		else
			% premultiply values
			wfa = w*len;
		end

		% add partial cell values to the output matrix
		for adx=1:nt2
		    % diagonal entry
		    I = 0.5*(wfa.*phi(:,adx))'*phi(:,adx);
		    nb = nb+1;
	   	    buf(nb,1) = T(idx,adx);
		    buf(nb,2) = T(idx,adx);
		    buf(nb,3) = I;
		    if (~flag)
			    % off diagonal entry
			    for bdx=adx+1:nt2
				I = (wfa.*phi(:,adx))'*phi(:,bdx);
				nb = nb+1;
				buf(nb,1) = T(idx,adx);
				buf(nb,2) = T(idx,bdx);
				buf(nb,3) = I;
			    end
		    end % not symmtric
		end % idx
	end
	% build matrix
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(P,1), size(P,1));
	% symmetric completition
	A = A+A';
end % assemble_1d_phi_phi

