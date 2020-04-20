% 2017-03-30 10:30:24.656509098 +0200
%% roots of linear functions
function r =roots1(t,y)
	if (nargin()<2)
		c  = t;
		r  = -c(:,2)./c(:,1);
	else
		r = t(:,1) - y(:,1).*(t(:,2)-t(:,1))./(y(:,2)-y(:,1));
	end
end

