% Tue Sep  6 09:28:07 CEST 2016
% Karl Kastner, Berlin
%
%% range of eigenvalues determined by the gershgorin circle theorem
function [lmin lmax ratio] = gershgorin_circle(A,pdflag)
	d  = diag(A);
	S  = sum(A-diag(D),2);
	lmin = d-abs(S);
	lmax = d+abs(S);
%	if (nargin() > pdflag)
%		lmin = 
%	end
	ratio = S./d;
end

