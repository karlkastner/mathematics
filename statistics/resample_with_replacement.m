function xx = resample_with_replacement(x)
	n = length(x);
	I = ones(n,1);
	N = (1:n)';
	for idx=1:n
		II = 1;
		for jdx=1:n
			if (idx~=jdx)
				II = kron(II,I);
			else
				II = kron(II,N);
			end
		end
		xx(:,idx) = x(II);
	end
end

