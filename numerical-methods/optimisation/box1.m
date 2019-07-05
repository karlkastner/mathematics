% Mon 15 Aug 17:03:04 CEST 2016
%% test objective function for optimisation routines 
function f=box1(c,x)
	f =   exp(-c(1)*x) ...
	    - exp(-c(2)*x) ...
	    - (exp(-x)-exp(-10*x));
end

