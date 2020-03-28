% Sun 15 Mar 14:06:14 +08 2020
% coefficients of the multiplicative interaction of two fourier series
% 0,1,2,3,4
function ab = fourier_multiplicative_interaction_coefficients(a,b,m,kmin)
	na  = size(a,2);
	nb  = size(b,2);
	if (nargin()<3)
		m = na+nb-1;
	end
	ab = zeros(size(a,1),m);
	if (issym(a))
		ab = sym(ab);
	end

	% the interactiom with mean (0 frequency) is special, as it yields the same frequency,
	% do not write it
	if (nargin()<4)
		kmin = 1;
	end

	for idx=kmin:na
	 for jdx=kmin:nb
		% plus
		if (idx+jdx-1<=m)
			ab(:,idx+jdx-1) = ab(:,idx+jdx-1)+a(:,idx).*b(:,jdx);
		end

		% minus
		if (jdx>idx)
			ab(:,jdx-idx+1) = ab(:,jdx-idx+1) + conj(a(:,idx)).*b(:,jdx);
		else
			ab(:,-jdx+idx+1) = ab(:,-jdx+idx+1) + (a(:,idx)).*conj(b(:,jdx));
		end
         end
	end
	ab=ab/2;
end

