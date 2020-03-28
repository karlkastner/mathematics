function [Esum Ebin] = entropyS(h,edges)

	n = size(h,1);
	% bin width
	% assume bin width of outer bins to be that of the interior bins
	w = edges(2:end)-edges(1:end-1);
	W = repmat(w(:)',n,1);

	% entropy of individual bins
	Ebin = -h.*log(h./W)/log(2);
	Ebin(isnan(Ebin)) = 0;
	% entropy of histograms
	Esum = sum(Ebin,2);
end
