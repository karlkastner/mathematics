% 2024-05-15 21:25:57.375524101 +0200
% Karl Kastner, Berlin
function [h,hh] = normal_mixture(x,p,mu,sd)
	p = p/sum(p);
	h = p(1)*normpdf(x,mu(1),sd(1));
	if (nargout>1)
		hh = zeros(length(x),length(p));
		hh(:,1) = h;
	end
	for idx=2:length(mu)
		hi = p(idx)*normpdf(x,mu(idx),sd(idx)); 		
	 	h  = h + hi;
		if (nargout>1)
			hh(:,idx) = hi;
		end
	end
end


