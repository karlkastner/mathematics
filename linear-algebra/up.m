% Thu 16 Nov 13:57:58 CET 2017
% Karl Kastner, Berlin
%% up-row by shifting rows down
%% linearly extrapolating or circular treatment of uppermost row
function x = up(x,iscircular)
	if (size(x,1)<2)
		warning('up on row vector');
	end
	if (nargin()<2 || ~iscircular)
	x = [2*x(1,:)-x(2,:);
              x(1:end-1,:)];
	else
	x = [ x(end,:);
              x(1:end-1,:)];
	end
end

