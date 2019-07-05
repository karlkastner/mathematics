% Mo 4. Jan 12:30:27 CET 2016
% Karl Kastner, Berlin
%% fir filter with some fancy extras
function [Yi, serr] = filterp1(X,Y,Xi,wfunc,efunc,nflag)
	if (nargin() <5)
		efunc = [];
	end
	if (nargin() <6)
		nflag = true;
	end

	Yi = []; %zeros(length(Xi),size(Y,2));
	serr = [];
	% for each output point
	for idx=1:length(Xi)
		% get input points in range
		% compute filter window
		% TODO this can be sped up by binary search
		w = rvec(wfunc(X,Xi(idx)));
		if (nflag)
			w = w/sum(w);
		end
%		[mv mdx] = max(w)
%		Xi(idx)
		
		% TODO GLS and linear order
		% filter mean
		fdx = find(w ~= 0);
		if (~isempty(efunc))
			Yi(idx,:) = efunc(w(fdx),Y(fdx,:));
		else
			% weighed mean
			Yi(idx,:) = w(fdx)*Y(fdx,:);
	%		% standard error
			res = ones(length(fdx),1)*Yi(idx,:) - Y(fdx,:);
			serr(idx,:) = sqrt(1/((sum(w).^2/sum(w.^2))-1)*w(fdx)*res.^2);
		end
	end
end

