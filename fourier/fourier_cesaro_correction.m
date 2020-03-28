% 2020-01-04 21:26:09.849753320 +0800
function c_ = fourier_cesaro_correction(c,nf)
	c_(1,:) = c(1,:);
	for k=1:(length(c)-1)/2
		c_(2*k,:)   = (1-(k-1)/(nf+1))*c(2*k,:);
		c_(2*k+1,:) = (1-(k-1)/(nf+1))*c(2*k+1,:);
	end
end

