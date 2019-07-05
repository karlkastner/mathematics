% 2016-11-22 16:11:13.487438288 +0100
%% mean of the left and right element
function x = lrmean(x)
	x = [(x(1,:));
             0.5*(x(1:end-2,:) + x(3:end,:));
	     (x(end,:));
            ];
	%x = 0.5*(x(1:end-1)+x(2:end));
end

