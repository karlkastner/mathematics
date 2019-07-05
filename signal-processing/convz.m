% 2015-04-10 11:18:58.736674352 +0200
function X = convz(X,f)
	X = cvec(X);
	f = cvec(f);
	X = [zeros(size(f)); X; zeros(size(f))];
	X = conv(X,f,'same');
	X = X(length(f)+1:end-length(f));
end

