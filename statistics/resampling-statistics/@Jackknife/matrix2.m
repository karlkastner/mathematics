% Fri Oct  4 21:16:59 UTC 2013
% Karl Kästner, Berlin

%% matrix of estimations for jacknive with two samples left out

function out = jackknife2(func, varargin)
	% get number of samples
	n = size(varargin{1},1);
	mdx = 1;
	% leave one sample out
	for idx=1:n-1
	 % leave another sample out
	 for jdx=idx+1:n
		% select subsets of the function arguments
		subset = [1:idx-1 idx+1:jdx-1 jdx+1:n];
		for kdx=1:length(varargin)
			v{kdx} = varargin{kdx}(subset);
		end % for kdx
		% run estimation with sample subset
		out_m = func(v{:});
		% preallocate memory at the first run
		if (1 == mdx)
			out = zeros(n*(n-1)/2,prod(size(out_m)));
		end
		% write output to the output array
		out(mdx,:) = out_m(:);
		mdx = mdx+1;
         end  % for jdx
	end % for idx
end % jackknife2

