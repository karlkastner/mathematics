% Fri 11 Oct 17:06:30 PST 2019
function test()


a = 2;
b = 3;
c = 5;
d = -1;
Xi = [0,8];
L = diff(Xi);
opt = struct();
opt.nx = 100;
%[x, y, out] = bvp2fdmc(@odefun, @(x,a,b) bcfun(x,a,b,Xi), Xi, opt);
[x, y, out] = bvp2c(@odefun, @(x,a,b) bcfun(x,a,b,Xi), Xi, opt);
[x_, y_, out] = bvp2fdm(@odefun, @(x,a,b) bcfun(x,a,b,Xi), Xi, opt);
y(:,2) = y_;

f = ode();
y(:,3) = f(L,a,b,c,d,x);

figure(1);
clf();
subplot(2,2,1)
plot(x,abs(y))
ylabel('|y|');
subplot(2,2,2)
plot(x,angle(y))
ylabel('arg(y)');

N=2.^(2:9);
res = [];
% convergence
for idx=1:length(N)
	opt.nx = N(idx)
	
	[x, y, out]   = bvp2c(@odefun, @(x,a,b) bcfun(x,a,b,Xi), Xi, opt);
	[x_, y_, out] = bvp2fdm(@odefun, @(x,a,b) bcfun(x,a,b,Xi), Xi, opt);
	y(:,2) = y_;
	y(:,3)        = f(L,a,b,c,d,x);

	res(idx,:) = rms(bsxfun(@minus,y(:,1:2),y(:,3)));
end

figure(2);
loglog(N,res)
legend('bvp2c','bvp2fdm');
ylabel('rms(y - y_{analytic}');

function c_ = odefun(x,y)
	if (nargin()<1)
		c_ = zeros(0,4);
	else
	 	c_ = ones(length(x),1)*[a,b,c,d];
	end
end

function [v,p,q] = bcfun(x,~,~,Xi)
	switch (x)
	case {Xi(1)}
		p = [1,0];
		q = [1,1];
		v = 1;
	case {Xi(2)}
		p = [1,0];
		q = [1,1];
		v = 0;
	otherwise
		error('here');
	end
end % bcfun
end

function f = ode()

syms y(x) a b c d x L ;
s=dsolve(a*diff(y,x,2) + b*diff(y,x,1) + c*y + d == 0, y(0)==1,y(L)==0) 
f = matlabFunction(s)
end

