% Sat  4 Jan 23:43:10 +08 2020
% leaving out n-samples
function [yi,si] = smooth_with_spline(y,d,wrap,imethod)
	if (isvector(y))
		y = cvec(y);
	end
	if (nargin<4)
		%imethod = 'spline';
		imethod = 'pchip';
	end
	n = size(y,1);
	k = ceil(n/d)+1;
	if (wrap)
		y  = [y;y;y];
	end
	yi = zeros(n,size(y,2));
	yi2 = zeros(n,size(y,2));
	for jdx=1:size(y,2)
	for idx=1:d
		id = idx + (-2:k)*d;
		if (~wrap)
			id(id<1) = [];
			id(id>n) = [];
			id_ = id;
			%id = max(1,id);
			%id = min(id,n);
			%id = unique(id);
		else
			id_ = n+id;
		end
		%yi_ =	interp1(id,y(id,jdx),(1:n)',imethod,'extrap');
		yi_ =	interp1(id,y(id_,jdx),(1:n)',imethod); %,'extrapval',NaN);
		%yi_(1:id(1)-1) = y(1:id(1)-1);
		yi(:,jdx) = yi(:,jdx) + yi_;
		yi2(:,jdx) = yi2(:,jdx) + yi_.*yi_;
	end
	end
	yi = yi/d;
	si = sqrt((yi2/d - yi.^2)/(d-1));
end

