% Thu Jul 11 15:10:16 UTC 2013
% Karl Kästner, Berlin

function [range, rdx] = select_range(val, n, thresh, dthresh, time)
	d_    = abs(conv(val,1/(2*n)*[-ones(n,1); ones(n,1)],'same')) < dthresh ...
		 & val > thresh;
	d     = diff(d_);
	range = [find(d == 1); find(d == -1)];
	rdx   = d_(range(:,1)+1);
end % select_range

