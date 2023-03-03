function y  = kurtnpdf(x,mu,s,beta)
	y = beta./(2*sigma.*gamma(1/beta))*exp(-abs((x-mu)./s).^beta);
end

