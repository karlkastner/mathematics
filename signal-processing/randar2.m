% 2015-08-17 17:44:42.414235480 +0200
% Karl Kastner, Berlin
%
%% generate ar2 process
function [x D s] = randar2(sigma,r1,r2,n,m)
	r = [r1 r2];
%	r_ = roots([1 -r1 -r2])
	s2 =  (1-r(2))/(1+r(2))*1/((1-r(2))^2 - r(1)^2);
	s = sqrt(s2);
	R = sigma*randn(n,m);
if (0)
	e = sqrt(1/s2)*R;
	x = zeros(n,m);
	x(1,:) = e(1,:);
	x(2,:) = r1*x(1,:) + e(2,:);
	for idx=3:n
		x(idx,:) = r1*x(idx-1,:) + r2*x(idx-2,:) + e(idx,:);
	end
else
	D = spdiags(ones(n,1)*[-r(2), -r(1), 1]*sqrt(s2),-2:0,n,n);
	x = D \ R;
%	D = D/sqrt(s2);
end
%	norm(x-X)
%	x = [x X];
end

