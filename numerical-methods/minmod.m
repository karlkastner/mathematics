function mm = minmod(a,b)
	mm = (a.*b > 0).*(   a.*(abs(a) < abs(b)) ...
                           + b.*(abs(b) <= abs(a)) );
end
%
%function mml = minmod_limiter(r)
%	mml = minmod(1,r); 
%end
%
%
%
%function vll = van_leer_limiter(r)
%	vll = (r + abs(r))./(1 + abs(r));
%end
