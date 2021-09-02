% Wed 18 Nov 11:52:58 +08 2020
% different batches are sometimes sieved on different sieves
% this 
%
%function hi = resample(h,e,ei)
function hi = resample(obj,ei)
	h  = obj.h;
	e  = rvec(obj.edge);
	ei = rvec(ei);
	% make sure that smallest sieve diameter is on the left hand side
	if (e(end)<e(1))
		h = fliplr(h);
		e = fliplr(e);
	end
	if (ei(end)<ei(1))
		flipi = true;
		ei = fliplr(ei);
	else
		flipi = false;
	end

	% cdf
	h = cumsum(h,2);

	% normalize
	c = [zeros(size(h,1),1), h./h(:,end)];
	
	% interpolate
	if (isvector(c))
		ci        = interp1(e,c',ei,'linear');
	else
		ci        = interp1(e,c',ei,'linear')';
	end

	fdx = ei < e(1);
	ci(:,fdx) = 0;
	fdx = ei > e(end);
	ci(:,fdx) = 1;
	

	% avoid loss of mass due to truncation at right end
	ci(:,end) = 1;
	
	hi = diff(ci,[],2);
	if (flipi)
		hi = fliplr(hi);
	end
end

