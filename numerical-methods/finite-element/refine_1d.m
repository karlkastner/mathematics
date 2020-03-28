% Mon Apr 30 12:12:57 MSK 2012
% Karl Kästner, Berlin

% P : set of points, including boundary
% T : elements (line segments)
% M : index of elements to be refined
function [P T] = refine_1d(P, T, M)
	lp = size(P,1);
	lt = size(T,1);
	lm = size(M,1);

	% allocate memory
	P = [P; zeros(lm,1)];
	T = [T; zeros(lm,2)];

	% split the selected elements
	for idx=1:lm
		a = T(M(idx),1);
		b = T(M(idx),2);
		lt = lt+1;
		lp = lp+1;
		P(lp,1) = 0.5*(P(a) + P(b));
		T(M(idx),:) = [a lp];
		T(    lt,:) = [lp b];
	end % for idx
end % refine_1d

