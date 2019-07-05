% 2015-06-28 13:04:36.220577551 +0200
% Karl Kastner, Berlin
%% cross correlation of two sampled ar1 processes
function [rho xcf xcf_] = xar1(res1,res2)
	if (nargin() < 2)
		res2 = res1;
	end
	n = size(res1,1);
	for idx=1:size(res1,2)
		x = xcorr(res1(:,idx),res2(:,idx),'coeff');
		xcf(:,idx) = x(n:end);
	end
	xcf_ = mean(xcf,2);
	fdx = find(abs(xcf_) < exp(-1),1);
	rho = exp((fdx-1) \ log(xcf_(fdx)));
	%rho = xcf_(1);
end

