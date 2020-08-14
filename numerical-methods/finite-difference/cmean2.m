% Sun 17 May 10:59:46 +08 2020
function x = cmean2(x, p, n, varargin)
	if (nargin()<2||isempty(p))
		p = 0.5;
	end
	if (nargin()<3 || isempty(n))
		n = 1;
	end
	for idx=1:n
	x = (  (1-p)*x ...
	     + p*0.5*(  cmean(x,1,1,1,varargin{:}) ...
	              + cmean(x,2,1,1,varargin{:}) ) ...
	    );
	end
end
