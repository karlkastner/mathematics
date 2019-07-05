% 2016-08-15 17:07:53.750865216 +0200
function f = box2(c,x)
	f = exp(-c(1)*x) - exp(-c(2)*x) - c(3)*(exp(-x) - exp(-10*x));
end
