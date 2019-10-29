% Sat Sep 24 17:55:19 MSD 2011
% Karl Kästner

% see knuth

% adds two doubles into one quad = s + e
function [s e] = add(a, b)
	s  = a + b;
	a_ = s - b;
	b_ = s - a_;
	a_err = a - a_;
	b_err = b - b_;
	e = a_err + b_err;
end % add

