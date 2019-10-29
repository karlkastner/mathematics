% kahan summation, pairwise summation

function [s e] = acc(a)
	s = 0;
	e = 0;
	for idx=1:length(a)
%		a_ = a(idx) - e;
%		s_ = s + a_;
%		e  = ((s_ - s) - a_);
%		s = s_;
		[s e_] = sum_(s, a(idx));
		e = e + e_;
		[s e] = sum_(s, e);
	end
	%e = -e;
end

