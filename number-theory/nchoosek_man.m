% 2014-12-11 14:49:24.795090057 +0100
% Karl Kastner, Berlin
%
%% vecotrised binomial coefficient
%% b = N!/K!(N-K)!
function B = nchoosek_man(N,K)
	% TODO computing the nominator completely can cause overflow
	%	factorial_man(N,N-K)./factorial(K);
	fdx = N >= K;
	B = zeros(size(fdx));
	B(fdx) = factorial(N(fdx))./(factorial(N(fdx)-K(fdx)).*factorial(K(fdx)));
end

