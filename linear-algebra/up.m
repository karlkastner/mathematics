% Thu 16 Nov 13:57:58 CET 2017
function x = up(x)
	if (size(x,1)<2)
		warning('up on row vector');
	end
	x = [2*x(1,:)-x(2,:);
              x(1:end-1,:)];
end

