function y = quadfilt1(x,n)
	if (isvector(x))
		x = cvec(x);
	end
	w = quadratwin(n);
	y = x;
	for idx=1:size(x,2)
	y(:,idx) = conv(x(:,idx),cvec(w),'same');
	end
end
