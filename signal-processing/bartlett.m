% Sun  8 May 16:41:35 CEST 2016
%% Effective sample size factor for bartlett window
%% c.f. thiebaux
%% c.f spectral analysis-jenkins, eq. (6.3.27)
%% c = acf
%% note: results seams always to be 1 tac too low
%% T : reduction factor for dof
%% for ar1 with a = rho^k = exp(-k/L), T = 2L
function [T, fxx] = bartlett(c)
	if (isvector(c))
		c = cvec(c);
	end
	m = length(c);
	fxx = 1./(2*pi)*(c(1) + 2*sum( (1 - (1:m-1)'/m ).*c(2:m)));
	% ratio of effective sample size
	% note, c is already normalised by sigma
	T = 2*pi*fxx;
end

