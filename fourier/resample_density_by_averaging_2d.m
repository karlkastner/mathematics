% 2025-03-11 14:22:36.530277810 +0100
function S = resample_density_by_averaging_2d(S,no)
	if (length(no)<2)
		no(2) = no(1);
	end
	% for first and second dimension
	for idx=1:2
		S = resample_density_by_averaging(S,no(idx),idx);
	end
end

