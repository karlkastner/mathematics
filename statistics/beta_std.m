% Wed 22 Sep 09:04:29 CEST 2021
function sd = beta_std(a,b)
	s2 = a*b./((a+b).^2.*(a+b+1));
	% asymptotic
	% lim b-> inf : a/b^2
	sd = sqrt(s2);
end

