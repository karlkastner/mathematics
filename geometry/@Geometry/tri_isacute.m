% Mo 2. Nov 18:11:55 CET 2015
% Karl Kastner, Berlin
%
%% flag acute triangles
% X : [x1,x2,x3] x-coordinates of triangle corners
% Y : [x1,x2,x3] y-coordinates of triangle corners
%
function [isacute Xc Yc] = tri_isacute(X,Y)
	[Xc Yc] = Geometry.tri_excircle(X,Y);
	%isacute = Geometry.inTriangle(X,Y,Xc,Yc);
	isacute = Geometry.inTriangle( [X(:,1),Y(:,1)] ...
					, [X(:,2),Y(:,2)] ...
					, [X(:,3),Y(:,3)] ...
					, [Xc(:,1),Yc(:,1)] ...
				     );
end

