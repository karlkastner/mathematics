% Thu 16 Apr 13:48:11 +08 2020
function e = entropy(p)
	p = p./sum(p);
	e = p.*log(p);
	e(0 == p) = 0;
	e = -sum(e);
end
