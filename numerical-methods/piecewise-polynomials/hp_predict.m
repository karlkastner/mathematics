% Sun Feb  1 17:31:40 CET 2015
% Karl Kastner, Berlin
%
%% predict with piecewise hermite polynomial
function val = hp_predict(x0,x1,c,x)
	n = 0.5*(length(c)-2);
	% segment length
	dx = (x1-x0)/n;
	% determine segment index of target points
	sdx = floor((1/dx)*(x-x0))+1;
	sdx = min(n,max(1,sdx));
	% shift and scale target point coordinates to unit interval
	x = (x - x0)/dx - sdx + 1;
	% vandermonde matrix of x
	X = vander_1d(x,3);
	% inverse of the Hermite matrix on the unit interval
        Hi = [1     0     0     0
              0     1     0     0
             -3    -2     3    -1
              2     1    -2     1];
	% convert vandermonde to Hermite matrix
	C = X*Hi;
	% predict
	val = zeros(size(x));
	for idx=1:length(val)
		val(idx) = C(idx,:)*c(sdx(idx)*2-1:sdx(idx)*2+2);
	end
	% TODO prediction error
end % hp_predict

