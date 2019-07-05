% Fri Oct  4 21:16:59 UTC 2013
% Karl Kästner, Berlin

%% matrix of estimation for leaving out two samples at a time
%
% TODO : this function could be implemented recursively to
%        implement the d-delete jackknife
% 
function M = matrix1_STATIC(func, varargin)
	% get number of samples
	n = size(varargin{1},1);
	% leave one sample out at a time
	M = [];
	mask = true(n,1);
	for idx=1:n
		% mask out current value
		mask(idx) = false;
		% select subsets of the function arguments
		for kdx=1:length(varargin)
			v{kdx} = varargin{kdx}(mask,:);
		end % for kdx
		% estimation with sample subset
		M(:,:,idx) = feval(func,v{:});
		% preallocate memory at the first run
		if (1 == idx)
			% the second dimension will be of length 1 if the
			% function returns scalars
			help = zeros([size(M) n]);
			help(:,:,1) = M;
			M = help;
		end % if
		% unmask current value
		mask(idx) = true;
	end % for idx
end % matrix1_STATIC

