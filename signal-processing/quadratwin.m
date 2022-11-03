% Tue  6 Sep 12:22:24 CEST 2022
function w = quadratwin(n)
	t = (0:n);
	t = t-mean(t);
	A = vander_1d(t,2);
	w = vander_1d(0,2)*inv(A'*A)*A';
	%plot(w)
end

