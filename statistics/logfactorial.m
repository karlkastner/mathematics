% 2018-12-10 16:51:27.183151982 +0100
% Karl Kastner, Berlin
%
%% approximate log of the factorial
%
function f = logfactorial(n)
	f = 1/2*log(2*pi*n) + n.*(log(n)-1);
end

