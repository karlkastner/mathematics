% 2011-09-05 16:58:42.000000000 +0200
% Karl Kastner, Berlin

function X = root4(P)
	a = P(:,1);
	b = P(:,2);
	c = P(:,3);
	d = P(:,4);
	e = P(:,5);

	alpha = -3/8*b.^2./a.^2 + c./a;
	beta  = b.^3./(8*a.^3) - b.*c./(2*a.^2) + d./a;
	gamma = -3*b.^4./(256*a.^4) + c.*b.^2./(16*a.^3) - b.*d./(4*a.^2) + e./a;

	% 8 y^3 + 8 pm^2 + (2p^2-8r)m - q^2 = 0
	p = -alpha.^2/12 - gamma;
	q = -alpha.^3/108 + alpha.*gamma/3 - beta.^2/8;

%	y = zeros(size(p));
	
%	if (0 == p)
	fdx = (0 == p);
	y(fdx) = -cbrt(q(fdx));
	y(~fdx) = NaN;
if (0)
%	else
		z_r     = -0.5*q; % ankathete
		z_i_sqr =  0.25*q.*q + (1.0/27.0).*p.*p.*p;
		gdx = (z_i_sqr >= 0) & fdx;
		%if (z_i_sqr >= 0)
%			s2 = -0.5*q(gdx) + sqrt( s1(gdx) );
%			u = cbrt( s2 );
%			y(fdx&gdx) = u - (1.0/3.0)*p(gdx)./u;
			y(fdx&gdx) = NaN;
		%else
		gdx = (z_i_sqr < 0) & fdx;
		r = sqrt(z_r(gdx).*z_r(gdx) - z_i_sqr(gdx)); % radius
			% equal results for + pi*2/3 and pi *4/3
			% cbrt(sqrt(-0.5*q + sqrt( s1 )) + cbrt(-0.5*q - sqrt( s1 )) )
			y(gdx) = cbrt(r).*2.0.*cos(1.0/3.0*acos(z_r(gdx)./r)); 
end	

	y = y - 5/6*alpha; % y is always real

	w = sqrt(alpha + 2*y);
	
	X(:,1) = -0.25*b./a + 0.5*( w + sqrt(-(alpha + 2*y) - 2*(alpha+beta./w)));
	X(:,2) = -0.25*b./a + 0.5*( w - sqrt(-(alpha + 2*y) - 2*(alpha+beta./w)));
	X(:,3) = -0.25*b./a + 0.5*(-w + sqrt(-(alpha + 2*y) - 2*(alpha-beta./w)));
	X(:,4) = -0.25*b./a + 0.5*(-w - sqrt(-(alpha + 2*y) - 2*(alpha-beta./w)));
end % qroot

