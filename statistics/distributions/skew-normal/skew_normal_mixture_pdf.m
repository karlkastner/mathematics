% 2024-05-15 21:25:57.375524101 +0200
% Karl Kastner, Berlin
function [h,hh] = skew_normal_mixture(x,p,mu,sd,sk)
	p = p/sum(p);
	h = p(1)*skewpdf(x,mu(1),sd(1),sk(1));
	if (nargout>1)
		hh = zeros(length(x),length(p));
		hh(:,1) = h;
	end
	for idx=2:length(mu)
		hi = p(idx)*skewpdf(x,mu(idx),sd(idx),sk(idx)); 		
	 	h = h + hi;
		if (nargout>1)
			hh(:,idx) = hi;
		end
	end
end


