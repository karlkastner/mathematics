function [Y Z] = adapt(X,c,e)
%e = 2;
d = X(end)/(e*(exp(e*X(end))-1)); X = d*e*sign(X).*(exp(e*abs(X)) - 1);

Z=X;
	Y(1) = X(1);
	for idx=2:length(X)
		f = exp(-Y(idx-1));
		h = c/sqrt(f);
		Y(idx) = Y(idx-1)+h;
	end
end

