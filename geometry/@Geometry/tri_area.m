% 2015-03-03 15:44:12.130265557 +0100
% Karl Kastner, Berlin
%
%% angle of a triangle
% function A = triarea(X,Y)

function A = triarea(X,Y)
	ax = X(:,2)-X(:,1);
	bx = X(:,3)-X(:,1);
	ay = Y(:,2)-Y(:,1);
	by = Y(:,3)-Y(:,1);

	A = 0.5*(ax.*by - ay.*bx);
end

