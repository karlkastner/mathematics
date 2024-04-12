% Sat  2 Mar 14:42:10 CET 2024
function x = cycle1(obj,b,x)
	obj.s(1).b = b;
	obj.s(1).x = x;
	obj.cycle(1);	
	x = obj.s(1).x;
end

