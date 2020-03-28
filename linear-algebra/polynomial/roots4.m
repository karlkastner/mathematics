% 2011-09-05 16:58:42.000000000 +0200
% Fri 27 Mar 10:48:15 +08 2020
% Karl Kastner, Berlin
%
% c.f. wikipedia Quartic_function#Solution
function r = roots4(c)
	r = roots4_(c(:,1),c(:,2),c(:,3),c(:,4),c(:,5));
end
function r = roots4_(a,b,c,d,e)	
	p = (8.*a.*c - 3.*b.^2)./(8.*a.^2);
	q = (b.^3 - 4.*a.*b.*c + 8.*a.^2.*d)./(8.*a.^3);
	d0 = c.^2 - 3.*b.*d + 12.*a.*e;
	d1 = 2.*c.^3 - 9.*b.*c.*d + 27.*b.^2.*e + 27.*a.*d.^2-72.*a.*c.*e;
	Q = cbrt(0.5.*(d1 + sqrt(d1.^2 - 4.*d0.^3)));
	S = 0.5.*sqrt(-2/3.*p + 1./(3.*a).*(Q + d0./Q));
	r = [(-b./(4*a) - S) + 0.5*sqrt(-4*S.^2 - 2*p + q./S)*[-1,1], ...
	     (-b./(4*a) + S) + 0.5*sqrt(-4*S.^2 - 2*p - q./S)*[-1,1] ];
end
