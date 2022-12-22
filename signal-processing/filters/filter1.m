% Do 11. Jun 10:13:11 CEST 2015
% Karl Kastner, Berlin
%% filter along one dimension
function [Xf S S1] = filter1(X,n,order)
	% TODO odd - even
	nx = size(X,1)
	xf = (-n/2:n/2)';
	Xf = NaN(size(X),class(X));
	S  = NaN(size(X),class(X));
	S1 = NaN(size(X),class(X));
	for idx=n/2+1:nx-n/2
		x = X(idx-n/2:idx+n/2,:);
		% regression matrix
		A = vander_1d(xf,order+1);
		% coefficients
		c = A(:,1:end-1)\x;
		% prediction
		Xf(idx,:) = c(1,:);
		xf_ = A(:,1:end-1)*c;
		% residual
		res = xf_ - x; %bsxfun(@minus,x,Xf(idx,:));
		% higher order contribution
		%res1 = xf*(1/(xf'*xf)*xf'*res);
		res1 = A(:,end)*(A(:,end)\res);
		% total standard error
		S(idx,:) = sqrt(mean(res.*res)/(n-order-1));
		% contribution by linear term
		S1(idx,:) = sqrt(mean(res1.*res1)/(n-order-2));
	end
end % function filter

	
