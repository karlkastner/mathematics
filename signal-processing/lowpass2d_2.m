
function x = lowpass2d_2(x,L)
	s = 5;
	if (L>0)
	xf = -s*L:s*L;
	f = exp(-1/2*(xf/L).^2)
	f = f/sum(f);
	x = conv2(x,f,'same');
	x = conv2(x,f','same');
	%x = x-conv2(x,f,'same');
	end
end
