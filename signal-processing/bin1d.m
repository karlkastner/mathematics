% Tue Jun 24 16:10:45 WIB 2014
% Karl Kastner, Berlin
%% bin values of v sampled at x into bins bounded by "edges"
%% apply function v to it
function vhist = bin1d(x,v,edges,func)
	if (nargin() < 4)
		func = @mean;
	end
	if (length(x) ~= length(v))
		error('length of x and v must match');
	end
	[n, xdx] = histc(x,edges);
	% TODO use accummarray
	for idx=1:length(edges)-1
		fdx          = find(xdx == idx);
		% note: mean returns consistently NaN for empty bins
		vhist(idx,:) = feval(func,v(fdx));
		%vhist(idx,:) = mean(v(fdx));
	end
end % bin1d

