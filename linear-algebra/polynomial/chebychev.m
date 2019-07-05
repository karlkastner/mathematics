% Fri  7 Apr 14:02:38 CEST 2017
% 8.15 in dronkers
%% chebycheff polynomials
function c = chebychev(x,n)
	if (issym(x))
		syms c
	else
		c=zeros(1,n);
	end
	c(1) = 1;
	if (n>0)
		c(2) = x;
	end
	for idx=2:n
		c(idx+1) = 2*x*c(idx) - c(idx-1);
	end
end

