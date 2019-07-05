% So 22. Nov 11:00:34 CET 2015
% Karl Kastner, Berlin
%
%% base point (fusspunkt), i.e. point on a line with shortest distance
%% to another point
% x0 : point to be projected on x1-x2
function [xb, s, ds] = base_point(x0,x1,x2)
	%xb(1,:) = x0(1)*dir(1)
	n = size(x0,2);
	d21 = bsxfun(@minus,x2,x1);
	d01 = bsxfun(@minus,x0,x1);
	s =    ( d21(1,:).*d01(1,:) + d21(2,:).*d01(2,:) ) ...
	    ./ (d21(1,:).^2 + d21(2,:).^2);
	xb      = zeros(2,n);
	xb(1,:) = x1(1,:) + s.*d21(1,:); 
	xb(2,:) = x1(2,:) + s.*d21(2,:); 
%	xb = x1 + bsxfun(@times,s,d21);

	% TODO s is not yet normalised
	if (nargout() > 2)
		dxy = x0-xb;
		ds  = hypot(dxy(1,:),dxy(2,:));
	end
end % base_point

