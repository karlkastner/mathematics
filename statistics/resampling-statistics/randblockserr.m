% Sat 14 May 14:34:39 CEST 2016
% Karl Kastner, Berlin
%
%% standard error of sequentilly correlated data by blocking
%% block length should be sufficiently larger than correlation length
%% and sufficiently smaller than data length
%% this uses a sliding block approach, which reduces the variation of the error estimate
%
function serr = blockserr(err,m)
	if (isvector(err))
		err = cvec(err);
	end
	serr = NaN(1,size(err,2));
	n  = size(err);
	n_ = round(n(1)/m);
	%% TODO this does not work, randomly picking samples does not reveal the correlation
	id = randi(n_,n(1),1);
	%id = sort(id);
	mu = [];
	% todo sliding
	% todo, better randomly permute instead of randomly choosing
	for idx=1:n_
		mu(idx,:) = mean(err(id==idx,:));
	end
	serr = std(mu/size(mu,1));
if (0)
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
end
	% reshape into blocks of length m
%	err = cvec(err);
%	err = reshape(err,m,[]);
	% compute block mean
%	err = mean(err);
	% compute error variance
end

