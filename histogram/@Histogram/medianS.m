% 2016-03-08 21:39:40.358145745 +0800
function me = medianS(h,edge)
	H      = Histogram.cdfS(h);
	centre = Histogram.centreS(edge);
	valid  = Histogram.validS(h);

	me = NaN(size(H,1),1);
	% TODO quick fix make cdf strictly monotoneous
	for idx=2:size(H,2)
		H(:,idx) = max(H(:,idx),H(:,idx-1)+sqrt(eps));
	end
	for idx=1:size(H,1)
		if (valid(idx))
			me(idx)  = interp1(H(idx,:)',centre,0.5,'linear',NaN)';	
		end
		% error
		% actually the error is proportional to the second derivate
		% nearest leads to uniform error betwenn 0 and binwidth, so average of l1 err is bw/2
		%me_(idx,1) = interp1(H(idx,:)',centre,0.5,'nearest',NaN)';
		%me_(idx,1) = interp1(H(idx,:)',centre,0.5,'cubic',NaN)';	
	end
	%err = me - me_;
	%me = interp1(centre(:),H',0.5,'linear',NaN)';	
end % medianS

