% Sat Jan 25 13:35:28 WIB 2014
% Karl Kastner, Berlin
%
%% project all points onto the cross section and assign them nz-coordinates
%%
%% transform coordinate into N-T reference
%% rotate coordinate, so that cross section goes along x-axis
%% then x and y are n and t respectively scaled by width
%% N and T coordinates
%alpha = atan2(obj.cc(2), 1);
%alpha = atan2(obj.cs.dir(2), obj.cs.dir(1));
%c = cos(alpha); s = sin(alpha);

function [N, T] = xy2nt(X,Y,xy0,dir,iflag);
	% normalize
	dir = dir/norm(dir);
	c = dir(1);
	s = dir(2);
	R = [ c, -s;
	      s,  c];

	if (nargin() < 5)
		iflag = false;
	end

	if (~iflag)
		R = R';
	end
	N = zeros(size(X));
	T = zeros(size(Y));
	for idx=1:size(X,2)
		%NT =  (R*[ X(:,idx)' - x0;
                %           Y(:,idx)' - y0])';
		if (~iflag)
			X(:,idx) = X(:,idx) - xy0(1);
			Y(:,idx) = Y(:,idx) - xy0(2);
		end
		NT = [X(:,idx),Y(:,idx)]*R';
%		NT =  ([ X(:,idx) - xy0(1),...
 %                        Y(:,idx) - xy0(2)]*R');
		N(:,idx) = NT(:,1);
		T(:,idx) = NT(:,2);
		if (iflag)
			N(:,idx) = N(:,idx) + xy0(1);
			T(:,idx) = T(:,idx) + xy0(2);
		end
	end
end

