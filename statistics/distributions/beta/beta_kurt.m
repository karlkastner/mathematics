function ku=beta_kurt(a,b)
	ku = 3 + 6*((b-a).^2.*(a+b+1) - a.*b.*(a+b+2))./(a.*b.*(a+b+2).*(a+b+3));
end
