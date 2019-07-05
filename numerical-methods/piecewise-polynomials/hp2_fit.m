% Sun Feb  1 16:45:56 CET 2015
% Karl Kastner, Berlin
%%
%% fit a hermite polynomial
%% coefficients are derivative free
%% x0  : left point of first segment
%% x1  : right point of last segment
%% n   : number of segments
%% x   : sample x-value
%% val : sample y-value
%% c   : coefficients (values at points, no derivatives)
function [c, serr] = hp2_fit(x0,x1,n,x,val)
	% there must be at least 3 intervals
	n = max(3,n);
	% filter
	fdx = isfinite(val);
	x   = x(fdx);
	val = val(fdx);
	% TODO this condition is not sufficient
	if (length(x)<n+1)
		c = NaN(n+1,1,class(val));
		serr = NaN(class(val));
		return;
	end

	% segment length
	dx = (x1-x0)/n;
	% determine segment for each source point
	sdx = floor((1/dx)*(x-x0))+1;
	sdx = min(n-1,max(2,sdx));
	% shift and scale source point coordinates to unit interval
	x = (x - x0)/dx - sdx + 1;

	% find and average identical points
	[x id ai] = unique(x);
	val = accumarray(ai,val) ./ accumarray(ai,ones(size(val)));
	% reorder sdx
	sdx = sdx(id);
	
	% vandermonde matrix of source points
	X = vander_1d(x,3);

	% inverse of the Hermite matrix on the unit interval
	% note: D cancels
	% D = diag([1 1 1/2 1/6]);
	% Hi = inv([  1 0 0 0;	   % f  left
        %               1 1 1 1;	   % f  right
        %               0 1 0 0;	   % f' left
        %               0 1 2 3]); % f' right
	Hi = [	     1     0     0     0
		     0     0     1     0
		    -3     3    -2    -1
		     2    -2     1     1 ];

	% note: dx became unity by transformation
	T = [ 0,1,0,0;
              0,0,1,0;
           -1/2,0,1/2,0;
           0,-1/2,0,1/2];

	% convert vandermonde to Hermite matrix
	C = X*Hi*T;

	% construct regression matrix
	% A sparse matrix could be used here, however, a full matrix is optimal
	% if the number of source points is larger than the number of segments,
	% which has to be the case to make the regression matrix non-singular
	A = zeros(length(val),n+1);
	for idx=1:length(val)
		A(idx,sdx(idx)-1:sdx(idx)+2) = C(idx,:);
	end
	% TODO use QR

	% only compute values at segments at points that are supported,
	% i.e. reduce A to non-singular sub-matrix
	valid = any(A)';
	c = NaN(n+1,1,class(val));
	% invert to find coefficients
	c(valid) = A(:,valid) \ val;
	% standard error
	res = A*c - val;
	serr = sqrt(res'*res/(length(val)-length(c)));
end % hp2_regress

