function pol()
clf
m=4;
n=2^10;
a=0.0; b=4.0;
X = linspace(a,b,n)';
Y = sqrt(sin(pi*X)); %exp(-X); %sin(0.5*pi*X);
Y = sqrt(X);
Y = (X.*exp(-X));
Y = sinh(X);
Y = log(X+1e-3);
Y = X.^3;
Y = 1./(X+1).^2;
Y = sin(X)./cos(X);
Y = exp(-X);
Y =sin(0.5*pi*X);
Y = (X.*exp(-X));
Y = 1./sqrt(X); %(1 + X.^2 + X.^3);
%Y = sqrt((X+1).*(X-1)); %(1+X).*(X-1).*sqrt(1 - (X-1).^2);
%Y = (1+X).*(sin(0.5*pi*X));
%Y = -(2*X/b).^(1/2).*(X-b);
Y = abs(2*X/b).^(1/1.5).*abs(X-b).^(1/1);
yb = 1.1*max(Y);

subplot(2,2,1)
plot(X,Y,'k','Linewidth',2);
 hold on
Z = lpol(X,Y,m);
plot(X,Z,'k--','Linewidth',2);
legend('u(x)', '\pi(x)');
set(gca,'xtick',0:m);
set(gca,'ytick',0:m)
xlabel('x')
%plot(X,Z);
norm(Z-Y)
axis equal; %axis square
ylim([0 yb]); xlim([a b])
print -deps ../img/fem-interpolation.eps
%axis equal;
% axis square

subplot(2,3,4)
plot(X,Y); hold on
Z = lpol(X,Y,2*m);
plot(X,Z);
norm(Z-Y)
axis equal; %axis square
%yylim([0 1.1]); xlim([a b])
%axis equal; % axis square

subplot(2,3,2)
plot(X,Y); hold on
Z = qpol(X,Y,1);
plot(X,Z);
norm(Z-Y)
axis equal; %axis square
%yylim([0 1.1]); xlim([a b])
%axis equal; %axis square

subplot(2,3,5)
plot(X,Y); hold on
Z = qpol(X,Y,2*m);
plot(X,Z);
norm(Z-Y)
axis equal; %axis square
%yylim([0 1.1]); xlim([a b])

f = @(c) norm(Y - cpol(X,Y,c));
fminbnd(f,1,n)

end

function Z = qpol(X,Y,k)
	n = length(X);
	for idx=1:k
		xl = round((idx-1)*n/k+1);
		xm = round((idx-0.5)*n/k);
		xr = round((idx)*n/k);
		C = [ 1 X(xl) X(xl).^2;
	      	      1 X(xm) X(xm).^2;
       		      1 X(xr) X(xr).^2 ] \ [Y(xl) Y(xm) Y(xr) ]';
		f = @(x) C(1) + C(2)*x + C(3)*x.^2; %-(x - X((idx-1)*n/k+1)).*(x - X(idx*n/k))
		Z(xl:xr,1) = f(X(xl:xr));
	end
end

function Z = lpol(X,Y,k)
	n = length(X);
	for idx=1:k
		xl = round((idx-1)*n/k+1);
		xr = round((idx)*n/k);
		C = [ 1 X(xl)
       		      1 X(xr) ] \ [Y(xl) Y(xr) ]';
		f = @(x) C(1) + C(2)*x; %-(x - X((idx-1)*n/k+1)).*(x - X(idx*n/k))
		Z(xl:xr,1) = f(X(xl:xr));
	end
end

function Z = cpol(X,Y,c)
	n = length(X);
	for idx=1:2
		if (1 == idx)
			xl = 1;
			xr = floor(c);
		else
			xl = floor(c)+1;
			xr = length(X);
		end
		C = [ 1 X(xl)
       		      1 X(xr) ] \ [Y(xl) Y(xr) ]';
		f = @(x) C(1) + C(2)*x; %-(x - X((idx-1)*n/k+1)).*(x - X(idx*n/k))
		Z(xl:xr,1) = f(X(xl:xr));
	end
end
