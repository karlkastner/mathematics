% Thu 16 Nov 13:57:58 CET 2017
% Karl Kastner, Berlin
% up-element, i.e. down shift of the matrix along column direction
% either circular or linearly extrapolation at boundary
function x = up(x,iscircular)
	if (size(x,1)<2)
		warning('down on row vector');
	end
	if (nargin()<2 || ~iscircular)
	x = [x(2:end,:);
	     2*x(end,:)-x(end-1,:)];
	else
		% TODO use circshift
	x = [x(2:end,:);
	     x(1,:)];
	end
end

