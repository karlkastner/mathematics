function p = generalized_logistic_pdf(x,mu,a,b)
%	p = -log(exp(b.*(a - mu)) + 1)./(exp(b.*(a - mu)) + 1).^x;
	%p = (a*b*exp(a*(mu - x)))./(exp(a*(mu - x)) + 1).^(b + 1);
	%c = 1-1./(1 + exp(a.*(x-mu))).^b;
	p=(a*b*exp(-a*(mu - x)))./(exp(-a*(mu - x)) + 1).^(b + 1);
end


