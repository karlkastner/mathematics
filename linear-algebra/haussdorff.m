% 2015-07-23 14:24:31.519238155 +0200

%% haussdorf dimension
%% box counting: count cectangles passed through by line (covered by polygon)
%%
%% Koch snow flake 3:4 -> 1.2619
%% Kantor set      2:3, (4:9) ->  0.6309
%% quadrat         4:2, 9:3, 16:4 -> 2
%%
% note: a line must always have a dimension < 2
% note that this definition is similar to sinuosity
%
% note powers a^n, b^n yield same result, but scales n*a, n*b not!
% (take care of overflow for large k and power)
function [n r s ar]=cantor(k,a,b)
	if (k<=0)
		n=1;
		r=1;
		s=1;
		ar=1;
	else
		[n r s ar]=cantor(k-1,a,b);
		n = a*n;
		r = b*r;
		s = s+n/r;
		ar = ar + n/(r*r);
%		s  = s + n/r;
%		n = 2*n;
%		r = 3*r;
	end
	% depth, haussdorff dimension, sum
	hd = log(n)/log(r)
	[s]
	[ar]
end

