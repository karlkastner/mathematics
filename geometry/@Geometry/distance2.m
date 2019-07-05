% 2014-11-18 17:57:26.947277912 +0100
% Karl Kastner, Berlin
%
%% euclidean distance between two points
%% this function requires a and be of equal dimensions, or the least the first pair or second pair to be a scalar
% function d2 = distance2(a,b);
function d2 = distance2(a,b); %,a2,b2)
	d2 = sum(bsxfun(@minus,a,b).^2,2);
%	if (nargin() == 2)
%		retval =  (a(1)-b(1)).^2 + (a(2)-b(2)).^2;
%	else
%		retval =  (a-a2).^2 + (b-b2).^2;
%	end
end % dist2

