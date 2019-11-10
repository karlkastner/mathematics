% 2015-07-13 15:10:40.175947031 +0200
% Karl Kastner, Berlin
%
%% hodges lehman location estimator
%%
%% Asymptotic rms efficency of location estimte:
%%      mean:          1 s/sqrt(n)
%%      hodges lehman: sqrt(pi/3)*s ~ 1.0233 s/sqrt(n)
%%      median:        pi/2 s/sqrt(n) ~ 1.25 s / sqrt(n)
function m = hodges_lehman_location(X,dim)
	if (nargin()<2)
		dim = 1;
		if (isvector(X))
			X = cvec(X);
		end
	end
	if (2 == dim)
		X = X';
	end

	m = zeros(1,size(X,2));
	for idx=1:size(X,2)
		XX = bsxfun(@plus,X(:,idx),X(:,idx)');
		fdx = triu(true(size(XX)));
		m(1,idx) = 0.5*nanmedian(XX(fdx));
		% nonzeros fails for zero values
		% nonzeros(triu(XX)));
	end
end % hodges_lehman

