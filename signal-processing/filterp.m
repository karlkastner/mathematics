% 2014-08-24 20:11:00.181678980 +0200
function yf = filterp(y,l,order)
	n = length(y);
	for idx=1:length(y)
		fdx = (max(1,idx-l):min(idx+l,n))';
		fdx_ = fdx - idx;
		A = [];
		for jdx=1:order+1
			A = [A fdx_.^(jdx-1)];
		end
		W = cos(0.5*pi*fdx_/(order+2));
		W = W/sum(W);
		W = diag(W);
		c = (A'*W*A)\(A'*W*y(fdx));
		yf(idx,1) = c(1);
	end
end

