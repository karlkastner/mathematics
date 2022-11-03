% 2021-09-29 15:16:20.323378274 +0200
function [a,b] = gam_moment2param(m,s)
	a = m.*m./(s.*s);
	b = s.*s./m;
end

