% 2016-11-22 16:11:13.487438288 +0100
% Karl Kastner, Berlin
%% single gaussian smoothing step with kernel 1/4*[1,2,1]
function x = cmean(x,keepends)
	if (nargin()<2)
		keepends = false;
	end
	if (~keepends)
	x = [0.5*(x(1,:)+x(2,:));
             0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
	     0.5*(x(end-1,:)+x(end,:));
            ];
	else
             x(2:end-1,:) = 0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
		%x = 0.5*(x(1:end-1)+x(2:end));
	end
end

