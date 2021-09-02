% Sun 15 Mar 14:06:14 +08 2020
% coefficients of the multiplicative interaction of two fourier series
% with integer-multiple frequencies k=0,1,2,...,n, in exponential form:
%
% f_a = Re( sum_k a_k*exp(i o_k t) ) 
%     = 1/2*sum( a_k exp(i k o t) + conj(a_k) exp(-i k o t))
%
% (a_k*exp(i k o t) + a_k' exp(-i k o t)  ))*(b_l*exp(i l o t) + conj(b_l) exp(-i o l t)
%   = a_k  b_l  exp(i ( k + l) o t)
%   + a_k  b_l' exp(i ( k - l) o t)
%   + a_k' b_l  exp(i (-k + l) o t)
%   + a_k' b_l' exp(i (-k - l) o t)
%
%   = a_k  b_l  exp(i ( k + l) o t)
%   + a_k  b_l' exp(i ( k - l) o t)
%   + (a_k b_l' exp(i (k - l) o t))'
%   + (a_k b_l  exp(i (k + l) o t))'

function ab = fourier_multiplicative_interaction_coefficients(a,b,m,kmin)
	na  = size(a,2);
	nb  = size(b,2);
	if (nargin()<3 || isempty(m))
		m = na+nb-1;
		%m = na+nb;
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
	%if (1==kmin)
	%end

	for idx=kmin:na
	 for jdx=kmin:nb
		% plus
		if (idx+jdx-1<=m)
			ab(:,+idx+jdx-1) = ab(:,+idx+jdx-1) + a(:,idx).*b(:,jdx);
		end

		% minus
		if (jdx>idx)
			ab(:,+jdx-idx+1) = ab(:,+jdx-idx+1) + conj(a(:,idx)).*b(:,jdx);
		else
			ab(:,-jdx+idx+1) = ab(:,-jdx+idx+1) + (a(:,idx)).*conj(b(:,jdx));
		end
         end
	end
	% as each term is halfed, the sum is quared, but each term comes as a pair
	ab=ab/2;
end

