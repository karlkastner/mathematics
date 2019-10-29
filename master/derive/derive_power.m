% derivation by matrix powers

%A = sqrt(D)*L*sqrt(D); A=A*A; A=0.5*(A+A'); z=sqrt(D)*A*(sqrt(D) \ y); plot(x(3:end-2),z(3:end-2))
%n=1000; x = 2*pi*sort((0:n-1)'/(n-1).*(1 + 1e-7*randn(n,1))); y = sin(x); [L D] = laplacian_non_uniform([2*x(1)-x(2); x; 2*x(end)-x(end-1)]); plot(x,[y -L*y L*(L*y)]); ylim([-1 1]), [norm(y+L*y)],[ norm(y-L*(L*y))]
%>> A = sqrt(D)*L*sqrt(D); A=A*A; A=0.5*(A+A'); z=sqrt(D)*A*(sqrt(D) \ y); plot(x(3:end-2),z(3:end-2));  ylim([-1 1])                                                                                                                                           
%>> n=1000; x = 2*pi*sort((0:n-1)'/(n-1).*(1 + 1e-5*randn(n,1))); x=((0:(n-1))'/n).^4; x=2*pi*x/max(x); y = sin(x); [L D] = laplacian_non_uniform([2*x(1)-x(2); x; 2*x(end)-x(end-1)]); plot(x,[y -L*y L*(L*y)]); ylim([-1 1]), [norm(y+L*y)],[ norm(y-L*(L*y))]


function [D4 A2_] = derive_power(x)
	if (nargin > 0 && ~isempty(x))
		xk=x(1);
		xl=x(2);
		xc=x(3);
		xr=x(4);
		xs=x(5);
		hk=xl-xk;
		hl=xc-xl;
		hr=xr-xc;
		hs=xs-xr;
	else
		syms xk xl xc xr xs
		syms hk hl hr hs
	end

	% lower order derivatives
	A1 = [ 1 0 0 0 0;
		-hl/(hk*(hk + hl)), -(hk - hl)/(hk*hl), hk/(hl*(hk + hl)) 0 0;
		0 -hr/(hl*(hl + hr)), -(hl - hr)/(hl*hr), hl/(hr*(hl + hr)) 0
		0 0 -hs/(hr*(hr + hs)), -(hr - hs)/(hr*hs), hr/(hs*(hr + hs));
		0 0 0 0 1];
 
	A2 = [ 1 0 0 0 0;
		2/(hk*(hk+hl))  -2/(hk*hl) 2/(hl*(hk+hl)) 0 0;
		0 2/(hl*(hl+hr))  -2/(hl*hr) 2/(hr*(hl+hr)) 0;
		0 0 2/(hr*(hr+hs))  -2/(hr*hs) 2/(hs*(hr+hs));
		0 0 0 0 1];

	% higher order derivatives by matrix powers
	A3 = A1*A2;
	A4 = A2*A2;

	D4 = [ A1(3,:)
	       A2(3,:)
	       A3(3,:)
	       A4(3,:) ];
end % function derive_power

