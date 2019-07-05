% Tue Jun 24 14:40:35 WIB 2014
% Karl Kastner, Berlin
%% bin values of V sampled at X and Y into the grid structured grid ex,ey
%% apply function func to all walues in the bin
%% func = mean : default
%% func = sum : non-normalized frequency histogram in 2D
function [B, ex, ey] = bin2d(X,Y,V,ex,ey,func)
	if (nargin() < 6)
		func = @mean;
	end
	n = 11;
	% automatic calculation of boundaries
	if (nargin() < 4)
		q = quantile(X,[0.05 0.95]);
		ex = linspace(q(1),q(2),n);
	end
	if (nargin() < 5)
		q = quantile(Y,[0.05 0.95]);
		ey = linspace(q(1),q(2),n);
	end

	% histogram in first dimension
	[void xdx] = histc(X,ex);
	% histogram in second dimension
	[void ydx] = histc(Y,ey);
	n=1;
	for idx=1:length(ex)
	 for jdx=1:length(ey)
		fdx = find(xdx == idx & ydx == jdx);
		B(idx,jdx) = func(V(fdx,:));
		%B(n,:) = mean(V(fdx,:),1);
		n=n+1;
	 end
	end
end

