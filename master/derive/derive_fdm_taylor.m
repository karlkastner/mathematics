% 2012-05-22 10:43:09 UTC
% Karl Kästner, Berlin

% derivatives finite difference coefficients by means of Taylor expansion
% order has to be in [2 3 4]
function [D D_] = derive_fdm_taylor(o,x)
	if (nargin > 1 && ~isempty(x))
		xk=x(1);
		xl=x(2);
		xc=x(3);
		xr=x(4);
		xs=x(5);

		hk=xl-xk;
		hl=xc-xl;
		hr=xr-xc;
		hs=xs-xr;
		hkl = hk+hl;
		hrs = hr+hs;
	else
		syms xk xl xc xr xs
		%syms hk hl hr hs hkl hrs
		hk=xl-xk;
		hl=xc-xl;
		hr=xr-xc;
		hs=xs-xr;
		hkl = hk+hl;
		hrs = hr+hs;
	end

	switch (o)
	case {2}
	A2 = [   hr 1/2*hr^2;
	       -hl 1/2*hl^2 ];
	B2 = [0 -1 1;
	      1 -1 0];
	S2=A2 \ B2;
	D = [ 0 S2(1,:) 0
	      0 S2(2,:) 0
	      zeros(2,5) ];
	case {3}
	A3 = [        -hl      1/2*hl^2   -1/6*hl^3       1/24*hl^4
	               hr      1/2*hr^2    1/6*hr^3       1/24*hr^4 
	          (hrs) 1/2*(hrs)^2    1/6*(hrs)^3  1/24*(hrs)^4 ];
	B3 = [ 1 -1 0 0;
               0 -1 1 0;
	       0 -1 0 1];
	D = A3 \ B3
	D = [zeros(5,1) [D; zeros(1,4)]];

	A3 = [   -(hkl) 1/2*(hkl)^2   -1/6*(hkl)^3  1/24*(hkl)^4
                  -hl      1/2*hl^2   -1/6*hl^3       1/24*hl^4
	               hr      1/2*hr^2    1/6*hr^3       1/24*hr^4 ]
	B3 = [ 1 0 -1 0 ;
 	       0 1 -1 0 ;
               0 0 -1 1 ];
	D_ = A3 \ B3
	D_ = [[D_; zeros(1,4)] zeros(5,1)];

	case {4}

	A4 = [   -(hkl) 1/2*(hkl)^2   -1/6*(hkl)^3  1/24*(hkl)^4
	              -hl      1/2*hl^2   -1/6*hl^3       1/24*hl^4
	               hr      1/2*hr^2    1/6*hr^3       1/24*hr^4
	          (hrs) 1/2*(hrs)^2    1/6*(hrs)^3  1/24*(hrs)^4 ];
	B4 = [ 1 0 -1 0 0;
	      0 1 -1 0 0;
	      0 0 -1 1 0;
	      0 0 -1 0 1];

	D = A4 \ B4;
	end % switch

end
%{
D42 	 = [
 (2*(xr*xs+(xs+xr)*xl+2*(3/2*xc-xs-xl-xr)*xc))               /((xc-xk)*(xl-xk)*(xk-xr)*(xk-xs))
 (2*(xr*xs+(xs+xr)*xk+2*(3/2*xc-xs-xk-xr)*xc))               /((xc-xl)*(xk-xl)*(xl-xr)*(xl-xs))
(2*(xr*xs+(xs+xr)*xl+(xs+xl+xr)*xk + 3*(2*xc-xr-xs-xk-xl)*xc))/((xc-xk)*(xc-xl)*(xc-xr)*(xc-xs))
 (2*(xl*xs+(xs+xl)*xk+2*(3/2*xc-xs-xk-xl)*xc))            /((xr-xc)*(xk-xr)*(xl-xr)*(xr-xs))
 (2*(xl*xr+(xl+xr)*xk+2*(3/2*xc-xr-xk-xl)*xc))               /((xc-xs)*(xk-xs)*(xl-xs)*(xr-xs))
]
%}
