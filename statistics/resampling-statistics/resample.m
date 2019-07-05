% 2015-08-22 15:59:57.584403410 +0200
%
%% resample a vector and apply function to it
%%
%% TODO, should be with replacement
%%
%% n  : number of samples
%% m  : number of subsamples
%% cx : maximum number of combinations
function ret = resample(fun,n,m,cx)
	% TODO, length of c can explode, better random choice if number of
	% combinations is large
	c = nchoosek((1:n)',m);
	% if the number of combination increases the maximum,
	% choose a random subsample
	if (length(c) > cx)
		p = randperm(length(c));
		c = c(p(1:cx));
	end
	ret = fun(c);
end

