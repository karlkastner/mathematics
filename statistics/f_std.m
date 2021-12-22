% Wed 22 Sep 09:05:28 CEST 2021
function sd = f_std(d1,d2)
	sd = sqrt(2*d2^2*(d1+d2-2)./(d1*(d2-2)^2*(d2-4)));
end

