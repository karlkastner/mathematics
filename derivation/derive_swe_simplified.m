%Wed  2 Oct 20:15:04 CEST 2024
syms s1 hl hc hr zl zc zr dx2 K ns1 hu hd zu zd s2 ns2
syms sl sr su sd nsr nsl nsd nsu bl br dx
h = hc;
z = zc;
		  syms br bl
		  al = K*(h+z-hl-zl)/dx;
		  ar = K*(hr+zr-h-z)/dx;
                  dz_dt = (   0.5*ar.*( (1-br).*h   + (1+br).*hr) ...
                            - 0.5*al.*( (1-bl).*hl  + (1+bl).*h) ... 
                          )/dx;

Jl = simplify(diff(dz_dt,hl))
Jc = simplify(diff(dz_dt,hc))
Jr = simplify(diff(dz_dt,hr))


if (0)
	hc = h;
	dzl = (h +z  - hl-zl);
	%sl  = dzl > 0;
	dzr = (hr+zr - h -z);
	%sr  = dzr > 0;

	dz_dt = (K/dx2)*(        sr.*hr.*dzr -   sl.*hc.*dzl ...
			    +    (nsr).*hc.*dzr - (nsl).*hl.*dzl ...
			   );
%	dz_dt = (K/dx2)*(    s1 .*( hr.*((hr+zr)-(hc+zc)) - hc.*((hc+zc)-(hl+zl))) ...
%			    + (ns1).*( hc.*((hr+zr)-(hc+zc)) - hl.*((hc+zc)-(hl+zl))) ...
%			   );


	dzd = (h +z  - hd-zd);
	%sd  = dzd > 0;
	dzu = (hu+zu - h -z);
	%su  = dzu > 0;

	dz_dt = 0.*dz_dt + (K/dx2)*(        su.*hu.*dzu -   sd.*hc.*dzd ...
			    +    (nsu).*hc.*dzu - (nsd).*hd.*dzd ...
			   );

%	dz_dt = 0.*dz_dt + (K/dx2)*(    s2 .*( hu.*((hu+zu)-(hc+zc)) - hc.*((hc+zc)-(hd+zd))) ...
%			    + (ns2).*( hc.*((hu+zu)-(hc+zc)) - hd.*((hc+zc)-(hd+zd))) ...
%			   );
Jd = simplify(diff(dz_dt,hd))
Jc = simplify(diff(dz_dt,hc))
Ju = simplify(diff(dz_dt,hu))

end

%	s2    = (hu > hd);
%	dz_dt = dz_dt + (K/dx2(2))*(    s2 .*( hu.*((hu+zu)-(hc+zc)) - hc.*((hc+zc)-(hd+zd))) ...
%			    + (~s2).*( hc.*((hu+zu)-(hc+zc)) - hd.*((hc+zc)-(hd+zd))) ...
%			   );


