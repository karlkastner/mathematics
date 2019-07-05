% Wed 14 Mar 13:53:41 CET 2018
%% padd values to x
% TODO make matrix save
function y=paddval(x,n,val)
	if (nargin() < 2)
		n = 1;
	end
	if (nargin() < 3)
		val = 0;
	end
	y                      = zeros(size(x)+2*n);
	y(n+1:end-n,n+1:end-n) = x;
	y(:,1:n) = val;
	y(:,end-n+1:end) = val;
	y(1:n,:) = val;
	y(:,end-n+1:end) = val;
end

