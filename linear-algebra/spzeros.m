% 2015-09-28 14:44:52.700634864 +0200
%% allocate a sparze matrix of zeros
function Z = spzeros(n,m,varargin)
	if (nargin()<2)
		m = n;
	end
	Z = spalloc(n,m,0); %,varargin{:});
end

