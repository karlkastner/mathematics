% 2024-05-24 16:51:30.683850792 +0200
function c = cauchycdf(x,mu,gamma)
	c = 0.5+1/pi*atan((x-mu)/gamma);
end


