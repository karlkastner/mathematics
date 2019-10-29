% 2012 May  8 20:22 (MSK)
% Karl Kästner, Berlin

function P = derive_nc_3d(n)


P = basis_2d(n);
P(:,2) = 1 - P(:,2);
plot(P(:,2),P(:,3),'*')
A = vander_2d(P(:,2:3),n)
C = inv(vpa(A))

syms x y
X = vander_2d([x y],n);
f = X*C

F = int(f,y)
%F = subs(F,y,1-x) - subs(F,y,0)
F = subs(F,y,x) - subs(F,y,0)
F = int(F,x)
F = subs(F,x,1 ) - subs(F,x,0)
F=F*factorial(2*(n-1))/2
factorial(2*n-2)/2
round(F)
%[num den] = rat(F)

sum(num./den)
sum(num./den .*(num < 0))
sum(num./den .*(num > 0))

A = [1 0 0; 1 0 1; 1 1 0];
P(:,2) = 1 - P(:,2);
P = P*inv(A)

end

% todo, mark interior and corner points
function P = basis_2d(n)
	m = 1;
	for idx=0:n
	 for jdx=0:n-idx
	   P(m,:) = 1/n*[n idx jdx];
	   m = m+1;
	 end
	end
end


%{
% 2d 3 point
A = sym([	1 0 0;
	      	1 1 0;
	      	1 0 1]')
C = inv(A)

syms x y
syms f1 f2 f3
f = [f1 f2 f3]*C*[1; x; y]
F = int(f,y)
F = subs(F,y,1-x) - subs(F,y,0)
F = int(F,x)
I = subs(F,x,1 ) - subs(F,x,0)

pause

A = sym([      	1 1/2 0;
		1 1/2 1/2
		1 0 1/2]')
C = inv(A)

syms x y
syms f1 f2 f3
f = [f1 f2 f3]*C*[1; x; y]
F = int(f,y)
F = subs(F,y,1-x) - subs(F,y,0)
F = int(F,x)
I = subs(F,x,1 ) - subs(F,x,0)

pause
 P = sym(0.5*[	2 0 0;
	      	2 2 0;
	      	2 0 2;
	      	2 1 0;
		2 1 1
		2 0 1]')
A = [P; prod(P(2:3,:)); P(2:3,:).^2]

C = inv(A)

syms f1 f2 f3 f4 f5 f6
syms x y

f = [f1 f2 f3 f4 f5 f6]*C*[1 x y x*y x^2 y^2].'

F = int(f,y)
F = subs(F,y,1-x) - subs(F,y,0)
F = int(F,x)
F = subs(F,x,1 ) - subs(F,x,0)

pause
addpath('../fem');

 P = sym(0.5*[	2 0 0;
	      	2 2 0;
	      	2 0 2;
	      	2 1 0;
		2 1 1
		2 0 1]')
A = [P; prod(P(2:3,:)); P(2:3,:).^2]

C = inv(A)

syms f1 f2 f3 f4 f5 f6
syms x y

f = [f1 f2 f3 f4 f5 f6]*C*[1 x y x*y x^2 y^2].'

F = int(f,y)
F = subs(F,y,1-x) - subs(F,y,0)
F = int(F,x)
F = subs(F,x,1 ) - subs(F,x,0)

pause
end
%}

%{
% derive (n-1)-th order accurate Newton-Cotes
% Quadrature rule on the triangle
function derive_nc_3d(n)

% poit coordinates
P = [ 1 0 0;
      1 1 0;
      1 0 0];
% triangle
T = [1 2 3];

% get support points of the nth-order Lagrangian basis functions
[P T] = promote_2d(P,T,n);

A = vander_2d(P(T(1,:),n);
C = inv(A);
F = integrate_2d(C)
w = (F*ones(n+1,1))'

% integrate into y direction
...
% if the triangle is [1 0 0; 1 0 1; 1 1 1],
% than x is y
% sum y terms to x terms and zero them
...
% this yields a 1-dimensional intergral only
% integrate into x-direction

end % derive_nc_3d()

function polyint_2d_x(C)
	n = length(C,2);
	[1 
	%[1 1 2 1 2 3 2 2 3
end
%}


