% 2017-11-11 13:19:48.594392312 +0100
%% lax wendroff scheme
%% second order accurate, but no tvd
%% this is effectively not a limiter
%% eq. 6.39 in randall, leveque
function phi = lax_wendroff(theta)
	phi = ones(size(theta));
end

