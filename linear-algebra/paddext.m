% Wed 14 Mar 13:53:41 CET 2018
%% padd values to vactor
%% not suitable for noisy data
%% order = 0 : constant extrapolation (hold)
%% order = 1 : linear extrapolation
% TODO make matrix save
function y=paddext(x,n,order)
	
%	if (nargin() < 2)
%		val = 0;
%	end
	if (nargin() < 2)
		n = 1;
	end
	if (nargin() < 3)
		order = 0;
	end
	y            = zeros(length(x)+2*n,1);
	y(n+1:end-n) = x;
	switch (order)
	case {0}
		y(1:n)         = y(n+1);
		y(end-n+1:end) = y(end-n-1);
	case {1}
		y(1:n)         = x(1)   - (n:-1:1)*(x(2)-x(1));
		y(end-n+1:end) = x(end) + (1:n)*(x(end)-x(end-1));
%	
%	if (isempty(val))
%		y(1,2:end-1) = x(1,:);
%		y(end,2:end-1) = x(end,:);
%		% this also fills the corner value
%		y(2:end-1,1)   = y(:,2);
%		y(2:end-1,end) = y(:,end-1);
	end
end

