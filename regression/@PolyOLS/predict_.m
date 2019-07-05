% 2015-06-22 08:58:21.966365114 +0200
function [Y A] = predict(obj,param,X);
	A = vander_1d(X,obj.order);
	Y = A*param;
end

