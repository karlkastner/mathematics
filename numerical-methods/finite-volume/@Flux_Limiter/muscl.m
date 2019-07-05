% Thu 16 Nov 09:28:45 CET 2017
%% muscl flux limiter
function phi = muscl(theta)
	phi = max(0,min(2,2*theta),1/2(1+theta));
end

