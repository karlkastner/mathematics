% Do 11. Jun 10:07:45 CEST 2015
%% polynomial filter,
%% can be achieved by iteratively processing the data with
%% a mean (zero-order) filter
function [X S sd] = polyfilt1(X,nf,order)
	nfi = nf/order;
	for idx=1:order-1
		X = meanfilt1(X,round(nfi));
	end
	[X S sd] = meanfilt1(X,round(nfi));
end

