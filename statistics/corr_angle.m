function r = corr_angle(a,b)
	sa = sin(a);
	ca = cos(a);
	mua = atan2(mean(sa),mean(ca));
	sb = sin(b);
	cb = cos(b);
	mub = atan2(mean(sb),mean(cb));

	sa_ = sin(a-mua);
	sb_ = sin(b-mub);
	
	r = mean(sa_.*sb_)./sqrt(mean(sa_.*sa_.*sb_.*sb_));
end
