% Sun 29 Oct 16:22:56 CET 2017
%% roots of quadratic function
%% c1 x^2 + c2 x + c3 = 0
function r = roots2(c)
	r_ = sqrt(c(:,2).^2 - 4*c(:,1).*c(:,3));
	r = [(-c(:,2) + r_)./(2*c(:,1)), ...
	     (-c(:,2) - r_)./(2*c(:,1))];
end


