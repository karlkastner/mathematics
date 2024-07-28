function p =  normalfoldedpdf(x,mu,sd)
	p = zeros(size(x));
	fdx = x>=0;
	p(fdx) = normpdf(x(fdx),mu,sd) + normpdf(x(fdx),-mu,sd);
end

