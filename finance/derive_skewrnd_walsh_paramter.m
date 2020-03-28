function derive_skewrnd_walsh_parameter()
	% note that the skewness before and after transformation of mean and variance remains the same,
	% this is by definition
	syms a sk;
	mu1 = sqrt(2/pi)*(1-2*a);
	mu2 = 1;
	s2  = mu2 - mu1^2;
	mu3 = 9007199254740992/5644425081792261 - (18014398509481984*a)/5644425081792261;
	%a = solve(sk  - (mu3 - 3*mu1*s2 - mu1^3)/s2^(3/2), a)
	%eq = simplify(sk*s2^(3/2)  - (mu3 - 3*mu1*s2 - mu1^3),'ignoreanalyticconstraints',true)
	eq = simplify(sk^2*s2^(3)  - (mu3 - 3*mu1*s2 - mu1^3)^2,'ignoreanalyticconstraints',true)
	eq = expand(eq)
	eq = collect(eq,'a')
	digits(16)
	vpa(eq)
	p = [
	(- 16.51278562979814*sk^2 - 66.05114251919258);
	(49.53835688939443*sk^2 + 198.1534275575777);
	(- 42.46927885241419*sk^2 - 221.7535614345337);
	(2.374629555837648*sk^2 + 113.2514102731045);
	(6.060327092646548*sk^2 - 25.56209068255566);
	(1.008750944333701*sk^2 + 1.961956805599696);
	0.04798261113971328*sk^2 - 0.04752993595256072
	];
	r = roots(p)

end
