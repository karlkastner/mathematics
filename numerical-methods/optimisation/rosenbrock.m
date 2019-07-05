% 2016-08-15 16:49:12.593353019 +0200
%% rosenbrock test function
function Q = rosenbrock(x)
	Q = 100*(x(2)^2 - x(2)).^2 + (1 - x(1)).^2
end
