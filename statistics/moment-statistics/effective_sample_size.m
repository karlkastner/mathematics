% 2016-05-08 12:32:47.275283001 +0200
% Karl Kastner, Berlin
%
%% effective sample size of the weighted mean of uncorrelated data
%% c.f. Kish
function neff = effective_sample_size(w)
	neff = sum(w).^2./sum(w.^2);
end

