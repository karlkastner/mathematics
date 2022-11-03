% 2021-06-23 18:56:45.535572579 +0200

function x = bandpass2d(x,p,k) %varargin)
	for idx=1:k
%	x = bandpass1d(x,p,1); %varargin{:});
%	x = bandpass1d(x.',p,1)'; %varargin{:}).';
	x = bp_diagonal(x,p);
	x = bp_diagonal2(x,p);
	end
end

function x = bp_diagonal2(x,p)
	x =   lp_diagonal2(x,p);
	x = x-lp_diagonal2(x,p);
end

function x = bp_diagonal(x,p)
	x =   lp_diagonal(x,p);
	x = x-lp_diagonal(x,p);
end

function x = lp_diagonal(x,p)
	p_ = sqrt(p);
	for idx=2:size(x,2);
		x(idx,2:end) = (1-p_)*x(idx,2:end) + p_*x(idx-1,1:end-1);
	end
end

function x = lp_diagonal2(x,p)
	p_ = sqrt(p);
	for idx=2:size(x,2);
		x(idx,1:end-1) = (1-p_)*x(idx,1:end-1) + p_*x(idx-1,2:end);
	end
end
