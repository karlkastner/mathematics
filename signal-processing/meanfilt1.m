% Wed Jan 28 17:38:56 CET 2015
% Karl Kastner, Berlin
%% moving average filter with special treatment of the boundaries
function [Y, S, sd] = meanfilt1(X,nf)
	if (~isscalar(nf))
		error('nf must be scalar');
	end
	if (isvector(X))
		X = cvec(X);
	end

	% treatment of boundary
	method = 0;

	n = size(X,1);
%	Y = zeros(size(X));
	S = zeros(size(X));
	% TODO, this is not central, should be 0.5 at both ends if order is even
	nf = min(nf,n); % was floor(nf/2)
	l = floor(nf/2)
	r = nf-l-1
	% round(nf/2)-1;
	% TODO this is asymmetric at the ends
	% implement alternative where the filter width is reduced at the ends

	% inner part
	% convolution filter
	if (1 == mod(nf,2))
		% odd
		window = ones(1,nf)/nf;
	else
		% even
		window = [0.5,ones(1,nf-1),0.5]/nf;
	end
%	for idx=l+1:n-r % was + 2
%		[Y(idx,:),S(idx,:)] = mean_man( ...
%			X(idx-l:idx+r,:)');
%		[idx-l,idx,idx+r]
%	end
	Y = conv2(X,window','same');

	for idx=1:l+0 % was 1
		switch (method)
		case {0} % decrease windown size at boundary,
			 % first point is not filtered
			% in this case the end is not filtered at all
			% TODO, this can overflow if nf>n/2
			[Y(idx,:),S(idx,:)] = mean_man( ...
				X(1:2*idx-1,:)');
			% [1,idx,2*idx-1]
		case {1} % linear prediction at boundary,
			 % filter is half as wide at boundary
			A = [ones(l+idx,1),(1:l+idx)'];
			c = A \ X(1:l+idx,:);
			Y(idx,:) = [1,idx]*c;
		case {2} % linear prediction at boundary
			% filter of same length
			A = [ones(nf,1), (1:nf)'];
			c = A \ X(1:nf,:);
			A = [1,idx];
			Y(idx,:) = A*c;
		end
	end

		switch (method)
		case {0} % reduce filter size (symmetric)
			for idx=r-1:-1:0
				Y(end-idx,:) = mean_man(X(end-2*idx:end,:)');
				% [n-2*idx,n-idx,n]
			end
		case {1}
		case {2}
			for idx=1:r
				A = [ones(nf,1),(1:nf)']
				c = A \ X(end-nf+1:end,:);
				A = [1,nf-r+idx]
				Y(end-r+idx,:) = A*c;
				n-r+idx
				
			end
		end

	% estimate the discretisation error
%	dY  = cdiff(Y);
	% sd2 = nf.^2/12*dY.^2;
%	sd = (12^-0.5)*nf*abs(dY);

%	f = ones(n,1)/n;
%	if (isvector(X))
%		X = conv(X,f,'same');
%	else
%		X_ = zeros(size(X));
%		nh = floor(n/2);
		% head and tail
%		for idx=1:nh
%			X_(idx,:)     = mean(X(1:n+idx,:));
%			X_(end-idx+1) = mean(X(end-n-idx+1:end,:)); 
%		end
		% centre
%		for idx=1:size(X,2)
%			X__(:,idx) = conv(X(:,idx),f,'same');
%		end
%		X_(nh+1:end-nh,:) = X__(nh+1:end-nh,:);
		% sum
		% sum of squares
		% mean
		% standard deviation
		% standard error
% filter      (rectangular, triangular, han ... )
% local fit   linar, quadratic
% 	      kriging
		% error estimation, this assumes normally distributed variables
		% sum of normal variables is again normally distributed with
		% mu = sum ai mui, s^2 = sum ai^2 si^2
		%s2 = (sum x^2 - n xbar^2) / (n-1)
%	end
end % meanfilt1

