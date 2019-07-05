% Sun Feb  1 17:31:40 CET 2015
% Karl Kastner, Berlin
%
%% prediction with pw hermite polynomial
%% c are values at support points
function [val, dval] = hp2_predict(x0,x1,c,x)
	n = (length(c)-1);
	% segment length
	dx = (x1-x0)/n;
	% determine segment index of target points
	sdx = floor((1/dx)*(x-x0))+1;
	sdx = min(n-1,max(2,sdx));
	% shift and scale target point coordinates to unit interval
	x = (x - x0)/dx - sdx + 1;
	% vandermonde matrix of x
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

	% predict
	val = zeros(size(x));
	for idx=1:length(val)
		val(idx) = C(idx,:)*c(sdx(idx)-1:sdx(idx)+2);	
	end
	% first derivative
	if (nargout()>1)
		dX = vanderd_1d(x,1,3);
		dC = dX*hi*t;
		for idx=1:length(val)
			dval(idx) = dC(idx,:)*c(sdx(idx)-1:sdx(idx)+2);	
		end
	end
	% TODO prediction error
end % hp_predict

