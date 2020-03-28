% Sun 12 Jan 19:40:33 +08 2020
function x     = skewrnd_walsh(mu,sd,sk,n)
	 if (length(n)<2)
		n = [n,1];
	 end
	 a     = parameter(sk);
	 x     = abs(randn(n));
	 p     = rand(n);
	 fdx   = p<a;
	 x(fdx)= -x(fdx);

	% mu = 2/sqrt(2*pi)*(-a + (1-a))
	 mu1 = sqrt(2/pi)*(1-2*a);
	 mu2 = 1;
	 % second moment is 1, independent of a
	 s2  = mu2 - mu1^2
	 s1  = sqrt(1-mu1^2);
	 %mu3 = 9007199254740992/5644425081792261 - (18014398509481984*a)/5644425081792261
	 %mu3/sd.^(3/2)

	 %var(x)
	 %sk  = (mu3 - 3*mu1*s2 - mu1^3)/s2^(3/2)
	 %sk  = (mu3 - 3*mu*sd^2 - mu^3)/sd^(3/2)

	 % TODO, skale for variance

	 % correct the mean
	 %skewness(x)
	 x = (x-mu1)*sd/s1+mu;
	 %skewness(x)
	 %pause
end



function a = parameter(sk)
	p = [
	(- 16.51278562979814*sk^2 - 66.05114251919258);
	(49.53835688939443*sk^2 + 198.1534275575777);
	(- 42.46927885241419*sk^2 - 221.7535614345337);
	(2.374629555837648*sk^2 + 113.2514102731045);
	(6.060327092646548*sk^2 - 25.56209068255566);
	(1.008750944333701*sk^2 + 1.961956805599696);
	0.04798261113971328*sk^2 - 0.04752993595256072
	];
	% there are only two real a, one leading to -sk and one to +sk
	a = roots(p);
	a = a(abs(imag(a))<sqrt(eps)*abs(a));
	%a = max(a); a(1)
	if (sk > 0)
		a = min(a);
	else
		a = max(a);
	end

%	mu1 = sqrt(2/pi)*(1-2*p3);
%	mu2 = 1;
%	sd^2 - 1*p2;
%	a =    40564819207303339161528614086293/81129638414606681695789005144064 ...
%	    - (5644425081792261*sd^3*sk)/18014398509481984 ...
%	    - (5644425081792261*mu^3)/18014398509481984 ...
%	    - (16933275245376783*mu*sd^2)/18014398509481984
%	mu3 = 9007199254740992/5644425081792261 - (18014398509481984*a)/5644425081792261;
%	s2  = sd^2;
%	mu1 = mu;
%	sk_   = (mu3 - 3*mu1*s2 - mu1^3)/s2^(3/2) 
end

