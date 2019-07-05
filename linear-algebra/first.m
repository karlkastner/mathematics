% 2017-07-23 22:54:35.626492549 +0200
% first n-elements of a vector
function x=first(x,n)
	if (nargin()<2)
		n=1;
	end
	x = x(1:n,:);
	
end
