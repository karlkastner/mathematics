% 2022-06-13 09:23:48.699939640 +0200
% TODO, different BCS
% for higher dimesions, values just sum up
function e = laplacian_eigenvalue(L,n,id)
	e = -4*((n-1)./L).^2.*sin(pi*id./(2*(n+1))).^2;
end

