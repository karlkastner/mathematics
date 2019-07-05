% Di 5. Jan 14:17:47 CET 2016
% Karl Kastner, Berlin
%% hanning window for change point detection
% TODO try bayesian approach
% for nonparametric break/change point detection

function w = hanchangewin(x,x0,range)
	w       = hanwin(x,x0,range);
	fdx     = (x < x0);
	w(fdx)  = -w(fdx)/sum(w(fdx));
	w(~fdx) =  w(~fdx)/sum(w(~fdx));
end

