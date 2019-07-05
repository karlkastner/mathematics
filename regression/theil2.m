% Tue  4 Jul 16:38:03 CEST 2017
%
%% Theil senn-estimator for two dimensions (glm)
%
% TODO, this should use gauss-seidel iteration
%       (fit iteratively all dimension by keeping the parameters of the other
%        dimensions constant)
function c = theil2(X,y,m)
	n = size(X);

	if (nargin()<3)
		m = n(1);
	end

	% subsample
	id = randi(n(1),[m,n(2)]);

	c = zeros(n(2),m);
	for idx=1:m
		A = X(id(idx,:),:);
		c(:,idx) = A\y(id(idx,:));
	end
	
	fdx = all(isfinite(c));
	
	c=nanmedian(c(:,fdx),2);
% -(x12*y2 - x22*y1)/(- x12^2 + x11*x22)
%  (x11*y2 - x12*y1)/(- x12^2 + x11*x22)
end

