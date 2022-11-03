% Mon 13 Jun 12:38:13 CEST 2022
% c.f. Processus Stochastiques et Mouvement Brownien
%
% note : this works and is exact in contrast to the other methods
%
function [y,C,sC] = brownian_noise_2d_cov(n,sy)
	if (length(n)<3)
		n(3)=1;
	end
	% generate the covariance matrix
	id = flat((0:n(1)-1)'*ones(1,n(2)));
	jd = flat(ones(n(1),1)*(0:n(2)-1));
	if (nargin()>2)
		jd = jd/sy;
	end
	C = zeros(n(1)*n(2));
	for idx=1:n(1)*n(2)
	 for jdx=1:n(1)*n(2)
		% note: scale factor of 0.5 applied below
		C(idx,jdx) = (  hypot(id(idx),jd(idx)) ...
                              + hypot(id(jdx),jd(jdx)) ...
                              - hypot(id(idx)-id(jdx),jd(idx)-jd(jdx)));
         end
        end
	%e = randn(n(1)*n(2),n(3))/pi;
	e = randn(n(1)*n(2),n(3));
	sC = sqrtm(C);
	y = sC*e;
	y = reshape(y,[n(1),n(2),n(3)]);
	% scaling
	y = (0.5/sqrt(n(1)))*y;
end % brownian_noise_2d_cov

