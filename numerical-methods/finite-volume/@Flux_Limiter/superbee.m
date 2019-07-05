% Fri 10 Nov 16:21:07 CET 2017
%% superbee limiter
function phi = superbee(theta)
	%phi = max(0,min(1,2*theta),min(2,theta));
	phi = max(max(0, min(1,2*theta)),min(2,theta));
end

