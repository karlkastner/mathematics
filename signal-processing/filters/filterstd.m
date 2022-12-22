% 2015-02-02 13:31:27.218375676 +0100
function x = filterstd(x,n)
	%fdx = isfinite(u);
	q = quantile(flat(x),[0.16 0.5 0.84]);
	%n = sum(fdx(:));
	%plot(sort(u(fdx)),(1:n)/n)
	%xlim([-2 2])
	%hold on
	%vertline(l(1))
	%vertline(l(2))
	%ylim([0 1]);
	s = 0.5*(q(3)-q(1));
	l = q(2)+n*s*[-1 1];
	mask = prod( double((x>l(1) & x < l(2))), ndims(x));
	for idx=1:size(x,ndims(x))
		x_ = x(:,:,idx);
		x_(~mask) = NaN;
		x(:,:,idx) = x_;
	end
end

