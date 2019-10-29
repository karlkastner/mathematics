% Fri Feb 24 18:40:52 MSK 2012
% Karl Kästner, Berlin

function A = assemble_3d_phi_phi(P, T, func, int)
	nt1 = size(T,1);
	nt2 = size(T,2);

	% degree of basis function polynomial (always an integer)
	% 1,4,10,20,35,... -> 0,1,2,3,4,...
	% lt2 = 1/6*nv^3 + nv^2 + 11/6*nv + 1
 	nv = round(1/(3*(3*nt2+(9*nt2^2-1/27)^(1/2))^(1/3))+(3*nt2+(9*nt2^2-1/27)^(1/2))^(1/3)-2);


	% get quadrature coordinates and weights
	[w b flag] = feval(int);
	if (0 ~= flag && isempty(func))
		% diagonal mass matrix
		mb0 = nt2;
	else
		% non diagonal mass matrix
		% number of buffer entries per triangle
		mb0 = nt2*(nt2+1)/2;
	end

	% exact number of elements in buffer after integration
	nb0 = nt1*mb0;

	% allocate memory
	buf = zeros(nb0,3);

	% integrate shares of each triangle
	nb = 0;
	for idx=1:nt1
		% element-matrix
		A = [1 P(T(idx,1),:);
	       	     1 P(T(idx,2),:);
		     1 P(T(idx,3),:);
		     1 P(T(idx,4),:)];

		% tetra volume
		volume = 1.0/6.0*abs(det(A));

		% quadrature points (centre-point)
		q = b*A;

		Va = vander_3d(P(T(idx,:),:), nv);
		% get test function coefficients
		C  = inv(Va);

		% evaluation of the hat-function at the quadrature point
		Vq  = vander_3d(q(:,2:4), nv);
		phi = Vq*C;

		% evaluate coefficient function at quadrature points
		% and premultipy values
		if (nargin() > 3 && ~isempty(func))
			wfa = volume*w.*feval(func, q(:,2:4));
		else
			wfa = volume*w;
		end


		% add share to each of the elements corner points
%		m = mb0*(idx-1);
		for adx=1:nt2
			% diagonal entry
		    	I = 0.5*(wfa.*phi(:,adx))'*phi(:,adx);
			nb = nb+1;
			buf(nb,1) = T(idx,adx);
			buf(nb,2) = T(idx,adx);
			buf(nb,3) = I;
			if (0 == flag || ~isempty(func))
				% off diagonal entries
				for bdx=adx+1:nt2
					%A(p(adx),p(bdx)) = A(p(adx),p(bdx)) + I;
			    		I = (wfa.*phi(:,adx))'*phi(:,bdx);
					nb = nb+1;
					buf(nb,1) = T(idx,adx);
					buf(nb,2) = T(idx,bdx);
					buf(nb,3) = I;
				end % for adx
			end
		end % for bdx
	end % for idx

	if (nb ~= nb0)
		warning('preallocation insufficient');
	end

	% build matrix
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(P,1), size(P,1));

	% symmetric completition
	A = A + A';
end % assemmble_3d_phi_phi

