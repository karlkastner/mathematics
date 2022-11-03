% 2015-06-28 11:14:13.608098441 +0200
% Karl Kastner, Berlin
%
%% autocorrelation of the residual
function [rho acf] = ar1(res)
	n = size(res,1);
	for idx=1:size(res,2)
		acf(:,idx) = autocorr(res(:,idx),n-1);
	end
	% TODO, what is better mean coefficient or coefficient of the mean?
	if (1)
	% coefficient of mean
	acf_ = mean(acf,2);
	fdx = find(acf_ < exp(-1),1);
	rho(1) = exp((fdx-1) \ log(acf_(fdx)));
	fdx = find(acf_/acf_(2) < exp(-1),1);
	rho(2) = exp((fdx-2) \ log((acf_(fdx)/acf_(2))));
	else
	% geomean of coefficient
		for idx=1:size(res,2)
			fdx = find(acf(:,idx) < exp(-1),1);
			rho(idx,1) = exp((fdx-1) \ log(acf(fdx,idx)));
			fdx = find(acf(:,idx)/acf(2,idx) < exp(-1),1);
			rho(idx,2) = exp((fdx-2) \ log(acf(fdx,idx)/acf(2,idx)));
		end
		rho
		% this favours small correlations
		rho = geomean(rho);
	end
end

