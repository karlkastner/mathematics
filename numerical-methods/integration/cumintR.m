% Wed  1 Nov 11:35:23 CET 2017
%% cumulative integral from right to left
%function inty = cumintR(y,dx);
function inty = cumintR(y,arg2);
	if (isvector(y))
		y    = cvec(y);
	end
	if (isscalar(arg2))
		% dx : arg2
		arg2 = -arg2;
	else
		% x : arg2 
		arg2  = cvec(arg2);
		%dx = cvec(dx);
	end
	inty = flipud(cumintL(flipud(y),flipud(arg2)));
end
