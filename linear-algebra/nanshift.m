% 2014-09-17 12:31:10.322336803 +0200
% Karl Kastner, Berlin
%
%% shift vector, but set out of range values to NaN
function x = nanshift(x,n)
	x = circshift(x,n);
	x(1:n) = NaN;
	x(end+n+1:end) = NaN;
end

