% 2017-03-30 10:30:24.656509098 +0200
%% roots of linear functions
function t0 =roots1(t,y)
	t0 = t(:,1) - y(:,1).*(t(:,2)-t(:,1))./(y(:,2)-y(:,1));
end

