% 2013-10-16 12:17:36.893154679 +0200
% Karl Kastner, Berlin
%
%% test if set B is subset of A in O(n)-runtime
%%
%% A : first set
%% B : second set
%% P : set of primes (auxiliary)
function r = issubset(A,B,P)
	p = prod(P(B));
	q = prod(P(A));
	r = (0 == mod(p,q));
end

