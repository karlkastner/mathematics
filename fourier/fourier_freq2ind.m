% 2023-02-02 10:06:11.533271532 +0100
function f = fourier_freq2ind(f,n)
if (0)
	if (mod(n,2) == 0)
	if (max(f(:)) >= n/2)
		max(f(:))
		error('f out of bounds')
	end
	if (min(f(:)) < -n/2)
		min(f(:))
		error('f out of bounds');
	end
	else
	if (max(f(:)) > (n-1)/2)
		max(f(:))
		error('f out of bounds')
	end
	if (min(f(:)) < -(n-1)/2)
		min(f(:))
		error('f out of bounds');
	end
	end
end
	fdx = f < 0;
	f(fdx)  = f(fdx) + n+1;	
	f(~fdx) = f(~fdx) + 1;
%	if (max(f(:))>n || min(f(:)<1))
%		erro('subscript out of range');
%	end
end

