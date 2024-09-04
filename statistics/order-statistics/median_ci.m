% Di 2. Feb 11:55:44 CET 2016
% Karl Kastner, Berlin
%
% method below is based on rank transformation:
% page 144, eq 21 in conover 1999 3rd ed 
% n.b. this is recommended for sample sizes larger 20
%
% further reading:
% hojo 1931, 1933
% chu 1955
% cadwell 1952
% david lane (hyperstat)
%
%% median and its confidence intervals under assumption of normality
%% se_me = sqrt(pi/2)*sd/sqrt(n) = 1.25331*sd/sqrt(n)
function [pl, pm, pu] = median_ci(n,P)
	%z = (norminv(1-P));
	z = (norminv(1-(1-P)/2));
	d = z*sqrt(n);

	if (1 == mod(n,2))
		rankl = 0.5*(n - d );
		ranku = 0.5*(n + d ); % was +1, but this leads to asymmetric (wrong) results
	else
		rankl = 0.5*(n - d ); % was -1, but that was to large for smal samples
		ranku = 0.5*(n + d ); % was +1
	end
	% limitation is necessary if there are insufficient many samples
	pl = max(0,rankl/n);
	pm = 0.5;
	pu = min(1,ranku/n);
end

