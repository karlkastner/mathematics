function df = double_factorial(n)
	df = zeros(size(n));
	for idx=1:numel(n)
		df(idx) = prod(n(idx):-2:1);
	end
end

