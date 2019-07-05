% Do 11. Feb 16:45:49 CET 2016
% Karl Kastner, Berlin
%
%% inverse edgeworth expansion
%% c.f. cornis fisher 1937
%% c.f. Rao 2010
%% c.f. 2.50 in hall
%% CHERNOZHUKOV 3.3
function q = edgeworth_quantile(mu,sigma,c3,c4,p,sk,ku,n)
%	p = norminv(p);
if (false)
	% This for some apparent reason does not work
	p1 = -1/6*(p.^2-1);
	dp1 = -1/3*p;
	p2 = -p.*(c4/24*(p.^2 - 3) + 1/72*c3^2*(p.^4 - 10*p.^2 + 15));
	p1_ = -p1;
	p2_ = p1.*dp1 - 1./2.*p.*p1.^2 - p2;	% was 1./p instead of 1.2p
	q = mu + sigma*(norminv(p) + p1_ + p2_);
	%q = mu + sigma*(p + p1_ + p2_);
else
	order = 3;

	% coris fisher 1960
	% Kolassa, p. 42
	% excess kurtosis
	ku = ku-3;

	k = ku/24;
	s = sk/6;
	ku = 1/6*(1 + 11*s^2 + sqrt(s^4 - 6*s^2+1));
	kl = 1/6*(1 + 11*s^2 - sqrt(s^4 - 6*s^2+1));
	if (k > ku)
		k = ku;
	elseif(k < kl)
		k = kl;
		%R3 = 0;
	end
	z  = norminv(p);
	R1 = z;
	R2 = sk/6*(z.^2-1);
	%R3 = (3*ku*(z.^3-3*z) - 2*sk^2*(2*z.^3 - 5*z))/72;
	if (order > 2)
		R3 = (1/24*ku*(z.^3-3*z) - 1/36*sk^2*(2*z.^3 - 5*z));
	else
		R3 = 0;
	end
	% TODO check if this in the region of validity, if not, skip kurtosis term
%	R2 = 0;
	q = R1 + R2/n.^0.5 + R3/n;
	q = mu+sigma/sqrt(n)*q;
end % if

end % function

