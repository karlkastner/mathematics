% Wed 15 May 10:48:16 CEST 2024
% Karl Kastner, Berlin
% moments of a bivariate gaussian mixture
function mu = g2_moments(par,n)
	p1 = par(1);
	p2 = 1-p1;
	mu1 = par(2);
	mu2 = par(3);
	s1 = par(4);
	s2 = par(5);
	
	mt = p1*mu1 + p2*mu2;
	m1 = mu1 - mt;
	m2 = mu2 - mt;

	% p1*m1 + p2*m2; % = 0
	% this is just the weighted sum of non-central gaussian moments with means shifted to m1 and m2
	mu = [mt,
             p1*(s1^2 + m1^2) + p2*(s2^2 + m2^2);
	     p1*m1*(3*s1^2 + m1^2) + p2*m2*(3*s2^2 + m2^2);
             (   p1*(3*s1*s1*(2*m1*m1 + s1*s1) + m1^4)  ...
	       + p2*(3*s2*s2*(2*m2*m2 + s2*s2) + m2^4) );
	     ( p1*(5*m1*s1*s1*(3*s1*s1 + 2*m1*m1)+m1^5) ...
	     + p2*(5*m2*s2*s2*(3*s2*s2 + 2*m2*m2)+m2^5) );
	     ( p1*(m1^6 + 15*m1^4*s1^2 + 45*m1^2*s1^4 + 15*s1^6) ...
	     + p2*(m2^6 + 15*m2^4*s2^2 + 45*m2^2*s2^4 + 15*s2^6) );
	     ];
	mu = mu(1:n);
end

