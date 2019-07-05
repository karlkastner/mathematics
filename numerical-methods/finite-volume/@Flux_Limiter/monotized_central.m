% 2017-11-11 13:20:44.224803094 +0100
%% monotonized central flux limiter
function phi = monotized_central(theta)
	phi = max(0, min(min(0.5+0.5*theta),2),2*theta);
end

