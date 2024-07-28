function p = logistic_pdf(x,mu,a)
%	p = -log(exp(b.*(a - mu)) + 1)./(exp(b.*(a - mu)) + 1).^x;
	%p = (a*b*exp(a*(mu - x)))./(exp(a*(mu - x)) + 1).^(b + 1);
	%p=(a*exp(-a*(mu - x)))./(exp(-a*(mu - x)) + 1);
	%p= a./(1 + exp(a*(mu - x)));
	%p=(a*exp(-a*(mu - x)))./(exp(-a*(mu - x)) + 1);
	p=(a*exp(-a*(mu - x)))./(exp(-a*(mu - x)) + 1).^2;
end

