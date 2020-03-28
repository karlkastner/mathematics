% 2014-11-19 14:44:05.782047545 +0100
% Karl Kastner, Berlin

% input h must be scaled by inverse bin width
function [maxdx mindx obj] = bimodes(obj)
	h = obj.h;
	n = length(h);
	[maxval maxdx] = max(h);
	mindx=[];
	minval_ = +Inf;
	% to find true side maxima this is 0
	maxval_ = 0;
	% right side
	for idx=maxdx+1:n
		if (h(idx) < minval_)
			minval_ = h(idx);
			mindx_  = idx;
		elseif(h(idx)-minval_ > maxval_)
			maxval_ = h(idx);
			maxdx(2) = idx;	
			% only accept minimum if it separates the maxima
			mindx    = mindx_;
		end
	end
	% reset minval, but not maxval
	minval_ = +Inf;
	% left side
	for idx=maxdx-1:-1:1
		if(h(idx)<minval_)
			minval_ = h(idx);
			mindx_  = idx;
		elseif(h(idx)-minval_ > maxval_)
			maxval_  = h(idx);
			mindx    = mindx_;
			maxdx(2) = idx;
		end
	end
end % bimodes

