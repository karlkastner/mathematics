% Mon Jun  4 22:27:34 MSK 2012
% Karl Kästner, Berlin

function y = mfunc_jdx(L, X, b)
	% (I - XX')L y = b
	y = L \ (b - (X'*(X*b)));
	y = L \ (b - (X'*(X*b)));
end

