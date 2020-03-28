% Fri 27 Mar 10:32:01 +08 2020
function c = roots2poly(r)
	n = size(r);
	if (n(2)<=1)
		c = [ones(n(1),1),-r];
	else
		c  = roots2poly(r(:,2:end));
		c  = [c,zeros(n(1),1)] + [zeros(n(1),1),-r(:,1).*c];
	end
end

