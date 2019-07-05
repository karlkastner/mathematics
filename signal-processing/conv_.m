% 2016-02-27 01:27:14.468176379 +0100
%% convolution of a with b
%
function c = conv_(a,b)
	na = length(a)
	nb = length(b);
	b  = flipud(b);
	c  = zeros(na-nb+1,1);
	for idx=1:na-nb+1
		c(idx) = sum(a(idx:idx+nb-1).*b);
	end
end
