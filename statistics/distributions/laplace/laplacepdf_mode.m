% 2024-06-30 11:06:48.285498931 +0200
function [fc,Sc] = laplacepdf_mode(f0,s)
	% syms s f f0 positive; S = laplacepdf(f,f0,s); fc=solve(diff(S,f),f), Sc = subs(S,f,f0)
	fc = f0;
	Sc = 1./(2*s);
end

