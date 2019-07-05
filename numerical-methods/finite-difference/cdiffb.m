% Thu Feb 26 11:23:04 CET 2015
% Karl Kastner, Berlin

%% differences of columns of X
%% degree  = 1 : central first order differences
%% degreee = 2 : central second order differences
%% TODO use difference matrix function for simplicity
function d = cdiff(X,degree,order)
	if (nargin() < 2 || isempty(degree))
		degree = 1;
	end
	if (nargin() < 3 || isempty(order))
		order = 2;
	end
	if (isvector(X) && isrow(X))
		X = X.';
		ir = true;
	else
		ir = false;
	end
	d = zeros(size(X));
	switch (degree)
	case {1}
		%d = [X(2) - X(1); 0.5*(X(3:end) - X(1:end-2)); X(end)-X(end-1)];
		if (order < 3)
			d(2:end-1,:) = 0.5*(X(3:end,:)-X(1:end-2,:));
		else % 5 point kernel
			d(3:end-2,:) = (1/12)*(X(1:end-4,:) - 8*X(2:end-3,:) + 8*X(4:end-1,:) - X(5:end,:));
			d(2,:)     = 0.5*(X(3,:)-X(1,:));
			d(end-1,:) = 0.5*(X(end-2,:)-X(end,:));
		end

		if (size(X,1) > 1)
			d(1,:)       = X(2,:)   - X(1,:);
			d(end,:)     = X(end,:) - X(end-1,:);
		end
	case {2}
		d(2:end-1,:) = (X(1:end-2,:)-2*X(2:end-1,:)+X(3:end,:));
		d(1,:)       = d(2,:);
		d(end,:)     = d(end-1,:);
	otherwise
		error('not yet implemented');
	end % switch
	if (ir)
		X = X';
	end % of else if isvector
end % cdiff

