% Mon 14 May 18:42:48 CEST 2018

% this will take at least half a day to complete

function derive_fourier_piecewise_constant()

%	for a constant outflow, phi must linearly increase along x or y
%	v = 0, x < x1, v = v0, phi = (x-x1) x1<x<x2, v=0, phi=1, x1 < x < x2

syms a b c d k y x p q r0 drdx drdy

%f = a*sin(k*x) + b*cos(k*x)
f    = [ (a*cos(k*x) + b*sin(k*x)), (c*cosh(k*y) + d*sinh(k*y)) ];
%fobj = (f - y - y1*x).^2
%fobj = (f - y0 - (y1-y0)/(x1-x0)*(x-x0)).^2
% fobj = (f - y0 - dy*x).^2
% TODO, linear
xy = [x,y];
yx = fliplr(xy);
coeff = [a;b;c;d]
for id=1:length(xy)
	coeff(id)
	fobj = (p*f(id) + q*diff(f(id),xy(id)) - r0 - drdx*x + drdy*y).^2;
	for kdx=1:2
	if (2==kdx)
		disp('for k=0')
		fobj=subs(fobj,k,0);
	end
	F = simplify(combine(expand(simplify(int(expand(fobj),xy(id)))),'sincos'))
	for jd=1:length(coeff)
		dF(jd,1) = diff(F,coeff(jd));
	end % for jd
	% TODO this has to be repeated for xy = 2*pi
	dF_ = subs(dF,yx(id),0)
	[A,rhs] = equationsToMatrix(dF,coeff);
	A = simplify(A)
	rhs = simplify(rhs)
	end % for kdx
	% coeff_  = solve(dF,coeff)
	%[A,rhs] = equationsToMatrix(dF,coeff)
	%F=simplify(combine(expand(simplify(int(expand(fobj),x))),'sincos'));
	%dF = [diff(F,a); diff(F,b)]; [A,rhs] = equationsToMatrix(dF,[a,b])
end % for id
%pause
%dF = [diff(F,a);
%      diff(F,b)]
%rhs = collect(collect(rhs,'cos(k*x)'),'sin(k*x)')
%simplify(A*k^2)
%simplify(rhs*k^2)
%pause
%F_ = @(x_,y_) subs(subs(F,x,x_),y,y_);

%Y = [0,y0,0]
%X = [x0,x1,x2,x3];
%FF = 0
%for idx=1:length(Y)
%	FF = FF + F_(X(idx+1),Y(idx)) ...
%	        - F_(X(idx),Y(idx));
%end
%d = [diff(FF,a);
%     diff(FF,b)];
%disp(d)
%s = solve(d,[a,b]);
%s.a = simplify(s.a);
%s.b = simplify(s.b);
%afun = matlabFunction(s.a);
%bfun = matlabFunction(s.b);

end % derive

