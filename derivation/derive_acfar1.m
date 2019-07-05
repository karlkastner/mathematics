% 2017-08-11 16:07:07.439778546 +0200

syms x1 x2 x3 x4 x5 x6 x7

x = [x1 x2 x3 x4 x5 x6 x7];
mu = mean(x)
eps = x-mu;

syms a
n = length(x);
for idx=1:n
	%a(idx,1) = expand(mean(eps(1:n-idx+1).*eps(idx:n)));
	%a(idx,1) = expand(1/n*sum(eps(1:n-idx+1).*eps(idx:n)));
	%a(idx,1) = expand(mean(x(1:n-idx+1).*x(idx:n)));
	a(idx,1)  = expand(mean(x(idx)*mu));
end
mu = expand(mean(x).^2);
syms r
for idx=1:length(a)
	a(idx)
	for i1=1:n
	 for i2=1:n
		a(idx) = eval(['subs(a(idx),x',num2str(i1),'*x',num2str(i2),',r^-',num2str(abs(i2-i1)),')']);
		if (1==idx)
		mu = eval(['subs(mu(idx),x',num2str(i1),'*x',num2str(i2),',r^-',num2str(abs(i2-i1)),')']);
		end
	end
	end
%	a(idx)
%	pause
end
a
mu
