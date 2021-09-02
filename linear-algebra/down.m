% Thu 16 Nov 13:57:58 CET 2017
function x = up(x)
	if (size(x,1)<2)
		warning('down on row vector');
	end
	x = [x(2:end,:);
	     2*x(end,:)-x(end-1,:)];
end

