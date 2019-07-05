% Thu Apr 28 04:52:01 MSD 2011
% So 7. Feb 13:25:02 CET 2016
% Karl Kästner
%
%% min-mod schock limiter
function phi = minmod(theta)
	phi = minmod_(1,theta);

function M = minmod_(A,B)
	M = (A.*B > 0).*((abs(A) < abs(B)).*A + (abs(A) > abs(B)).*B);
end
end % minmod


