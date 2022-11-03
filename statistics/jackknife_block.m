% 2021-09-10 19:34:56.107905873 +0200
% c.f. busing 1999
function [mu,sd] = jackkife_block(fun,y,m)
	n  = length(y);
	np = floor(n/m); 
	% estimate with all samples
	ht_n    = fun(y);
	% "delete" (window) samples and estimate
	for jdx=1:m
		y_ = y;
		y_((jdx-1)*np+(1:np)) = 0;
		ht(:,jdx) = fun(y_);
	end
	% average
	bt_m = mean(ht,2);
	% jacknive estimate
	mu = m*ht_n - (m-1)*bt_m;
	% pseudo values
	% tt(:,j) = g*ht_n - (g-1)*ht(:,j);
	% variance of the estimate
	s2 = (m-1)*mean((ht - bt_m).^2,2);
	sd = sqrt(s2);
end


