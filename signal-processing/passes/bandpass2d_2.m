function x = bandpass2d_2(x,L)
	%x = conv2(x',f,'same')';
	%x = x-conv2(x',f,'same')';
	x = lowpass2d_2(x,L);
	x = x-lowpass2d_2(x,L);
end

