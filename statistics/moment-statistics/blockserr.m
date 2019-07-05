% Sat 14 May 14:34:39 CEST 2016
% Karl Kastner, Berlin
%
%% estimate the standard error of potetially sequentilly correlated data
%% by blocking
%% block length should be sufficiently larger than correlation length
%% and sufficiently smaller than data length
%% this uses a sliding block approach, which reduces the variation of the error estimate
function serr = blockserr(err,m)
	if (isvector(err))
		err = cvec(err);
	end
	serr = NaN(1,size(err,2));
	% TODO circulate
	for idx=1:size(err,2)
		% compute sliging mean of blocks of length m
		win  = ones(1,m)/m;
		% average over blocks
		errc = conv(err(:,idx),win,'valid');
		% compute std
		sd = std(errc);
		% number of samples after blocking (can be a fraction)
		ns = length(errc)/m;
		%ns = size(err,1)/m;
		% standard error
		%serr(idx) = sd/sqrt(ns);
		serr(idx) = sd/sqrt(ns-1);
	end
	% reshape into blocks of length m
%	err = cvec(err);
%	err = reshape(err,m,[]);
	% compute block mean
%	err = mean(err);
	% compute error variance
end

