% 2017-03-14 15:36:24.856216337 +0100
%% make x orthogonal to Y
function x = orthogonalise(x,Y,flag)
	Y = orth(Y);
	x = x-Y*(Y'*x);
	if (nargin()>2&&flag)
		x=x/norm(x);
	end
end

