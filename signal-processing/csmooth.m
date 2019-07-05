% 2017-05-08 13:26:24.568005932 +0200
% note: smoothing radius ~ sqrt(n)
%% smooth recursively with [1,2,1]/4 kernel
function x = csmooth(x,n,p,circ)
	if (nargin()<2 || isempty(n))
		n = 1;
	end
	if (nargin() < 3 || isempty(p))
		p = 1;
	end
	if (nargin() < 4 || isempty(circ))
		circ = false;
	end
	for idx=1:n
		if (~circ)
		x = (1-p)*x ...
                     + p*[0.5*(x(1,:) + x(2,:));
	                  0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
                          0.5*(x(end-1,:)+x(end,:))];
		else
		x = (1-p)*x ...
                     + p*[0.25*x(end,:) + 0.5*x(1,:) + 0.25*x(2,:);
	                  0.25*(x(1:end-2,:) + 2*x(2:end-1,:) + x(3:end,:));
                          0.25*x(end-1,:)+0.5*x(end,:)+0.25*x(1,:)];
		end
	end % for idx
end % csmooth

