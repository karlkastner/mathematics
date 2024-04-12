% Wed 14 Feb 14:40:00 CET 2024
% c.f Johnson 1989
% note: this is of little use when the matrix is indefinite, as the lower bound
% is usually negative in this case and thus does not give a non-trivial bound
% between zero and the smallest eigenvalue
function sv_min = estimate_sv(A)
	aA = abs(A);
	aD = diag(aA);
	R  = A-diag(diag(A));
	aR = abs(R);
	sv_min = min(2*diag(aA) - 0.5*(sum(aA,2)+sum(aA)'))
	sv_min = min(diag(aA) - 0.5*(sum(aR,2)+sum(aR)'))
	
	r = sum(aR,2);
	c = sum(aR,1)';
	sv_min = min(aD - max(r,c))
end

