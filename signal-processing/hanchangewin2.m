% Di 5. Jan 15:22:31 CET 2016
% Karl Kastner, Berlin
%% nanning window for chage point detection
% for nonparametric break/change point detection
function w = hanchangewin(x,x0,range)
	%w = sin(2*pi*(x-x0)/range).^2 .* sign(x-x0) .* (abs(x-x0)<range);
	w = sin(2*pi*(x-x0)/range).^2 .* (abs(x-x0)<range);
%	[x sdx] = sort(x);
%	sdx(sdx)=(1:length(x))';
%	w1 = hanwin(x,x0,range);
%	fdx = x < x0;
%	w2 = -(fdx) + (~fdx);
%	w = conv(w1,w2,'same');
	fdx = (x > x0);
	w(fdx)  =-w(fdx)/sum(w(fdx));
	w(~fdx) = w(~fdx)/sum(w(~fdx));
%	clf
%	[x sdx] = sort(x);
%	plot(x,w)
%	vline(x0)
%	pause
%	w = w(sdx);
end

