% 2014-12-18 10:30:02.050309598 +0100
%
%% limit a by lower and upper bound
function a = limit(a,l,u)
	a = min(max(a,l),u);
end

