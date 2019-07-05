% Wed 14 Mar 13:53:41 CET 2018
%% padd values at end of x
% TODO make matrix save
function y=paddval1(x,n,val)
	if (nargin() < 2)
		n = 1;
	end
	if (nargin() < 3)
		val = 0;
	end
	y                      = zeros(length(x)+2*n,1);
	y(n+1:end-n) = x;
	y(1:n) = val;
	y(end-n+1:end) = val;
end

