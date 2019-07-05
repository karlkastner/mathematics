% 2016-03-04 12:36:46.874851789 +0100
% Karl Kastner, Berlin

%% pdf from sample distribution
%% Note: better use kernal density estimates
function pdf = pdfsample(x,xs)
	xs = sort(xs);
	n = length(xs);
	cdfs = (1:n)/(n+1);
	cdf = interp1(xs,cdfs,x,'linear');
	pdf = cdiff(cdf)/((x(2)-x(1))); % why not 2h *
	pdf(~isfinite(pdf)) = 0;
end

