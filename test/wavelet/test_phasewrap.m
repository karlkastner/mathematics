% 2016-07-11 11:13:36.791266180 +0200

	fdx = phi < -pi;
	phi(fdx) = 2*pi + phi(fdx);
	fdx = phi > pi;
	phi(fdx) = phi(fdx)-2*pi;

