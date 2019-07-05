% Fri Feb 20 11:22:20 CET 2015
% Karl Kasntner, Berlin
%% predict values and confidence intervals with the Theil-Sen method
%
function [Y S L U obj] = predict(obj,X)
	if (isvector(X) && isrow(X))
		X=cvec(X);
	end
	one = ones(size(X));
	Y = one*obj.param(1,:) + X*obj.param(2,:);
	% this is a rather pessimistic confidence interval,
	% as it yields the worst case combination of slope and intercept
	l1 = (X-obj.qx(3))*obj.qs(3,:);
	l2 = (X-obj.qx(1))*obj.qs(1,:);
	L  = one*obj.qy(1,:) + min(l1,l2);
	u1 = (X-obj.qx(1))*obj.qs(3,:);
	u2 = (X-obj.qx(3))*obj.qs(1,:);
	U = one*obj.qy(3,:) + max(u1,u2);
	S = U-L;
end % predict

