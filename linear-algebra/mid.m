% 2015-04-07 18:35:59.479232653 +0200
% 2016-09-02 23:03:31.679656300 +0200
% Mon 29 Jan 14:31:23 CET 2018
% 2018-03-22 10:02:24.961332914 +0100
%% mid point between neighbouring vector elements
function x = mid(x,dim)
	if (isvector(x))
		x = 0.5*(x(1:end-1)+x(2:end));
	else
		if (nargin()<2)
			dim = 1;
		end
		switch (ndims(x))
		case {2}
			if (1 == dim)
				x = 0.5*(x(1:end-1,:)+x(2:end,:));
			else
				x = 0.5*(x(:,1:end-1)+x(:,2:end));
			end
		case {3}
			switch (dim)
			case {1}
				x = 0.5*(x(1:end-1,:,:)+x(2:end,:,:));
			case {2}
				x = 0.5*(x(:,1:end-1,:)+x(:,2:end,:));
			case {3}
				x = 0.5*(x(:,:,1:end-1)+x(:,:,2:end));
			end
		case {4}
			switch (dim)
			case {1}
				x = 0.5*(x(1:end-1,:,:,:)+x(2:end,:,:,:));
			case {2}
				x = 0.5*(x(:,1:end-1,:,:)+x(:,2:end,:,:));
			case {3}
				x = 0.5*(x(:,:,1:end-1,:)+x(:,:,2:end,:));
			case {4}
				x = 0.5*(x(:,:,:,1:end-1)+x(:,:,:,2:end));
			end
			
		end
	end
end
