% 2016-03-08 15:04:46.362525316 +0100
function q = quantileS(h,edge,p)
	H      = Histogram.cdfS(h);
	centre = Histogram.centreS(edge);
	% append 1
	H(:,end+1) = 1;
	centre(end+1) = 2*centre(end)-centre(end-1);

	q = NaN(size(H,1),length(p));

	% make cdf strictly monotoneous
	for idx=2:size(H,2)
		H(:,idx) = max(H(:,idx),H(:,idx-1)+sqrt(eps));
	end
	valid = Histogram.validS(h);
	for idx=1:size(H,1)
		if (valid(idx))
			% compute quantile
			q(idx,:) = interp1(H(idx,:),centre,p,'linear',NaN);
		end
	end
end
