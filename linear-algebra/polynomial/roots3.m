function r = roots3(c)
	d0 = c(:,2).^2 - 3*c(:,1).*c(:,3);
	d1 = 2*c(:,2).^3 - 9*c(:,1).*c(:,2).*c(:,3) + 27*c(:,1).^2.*c(:,4);
	% TODO sign choice inside sqrt so that C is not zero
	C  = cbrt(0.5*(d1 + sqrt(d1.^2 - 4*d0.^3)));
	e = (0.5*(-1 + sqrt(-3))).^(0:2);
	r = -1./(3*c(:,1)).*(c(:,2) + C*e + d0./(C*e));
end

