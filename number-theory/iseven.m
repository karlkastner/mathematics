% 2015-06-04 16:26:35.738518489 +0200
%
%% true for even numbers in X
%
function X = iseven(X)
	X = (mod(X,2) == 0);
end
