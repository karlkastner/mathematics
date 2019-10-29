% Wed May  2 14:17:13 MSK 2012
% Karl Kästner, Berlin

% TODO use int array with 2 or 3 elements, as this scheme overflows in 3d

function key = hashkey(a,b,c)
%	base = 8388607;
%	if (a > base || b > base)
%		error('haskey','here');
%	end
%	if (a < b)
%		key = a*base + b;	
%	else
%		key = b*base + a;
%	end


	if (nargin > 2)
		A = sort([a b c]);
		key = A(1)*2^36 + A(2)*2^18 + A(3);
	else
		A = sort([a b 0]);
		key = A(1)*2^36 + A(2)*2^18 + A(3);
	end

%	if (nargin > 2)
%		key = int32([a b c]');
%	else
%		key = int32([0 a b]');
%	end
%	key = sort(key);
end

