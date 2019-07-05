% So 31. Mai 12:41:54 CEST 2015
% Karl Kastner, Berlin
%
%% normalise a vector or the columns of a matrix
%% note that the columns are independently normalised, and hence not necessarily
%% orthogonal to each other use the gram schmidt algorithm for this (qr or orth)
function x = normalise(x)
	if (isvector(x))
		% normalise vector
		x=x/norm(x);
	else
		% normalise columns
		for idx=1:size(x,2)
			x(:,idx) = x(:,idx)./norm(x(:,idx));
		end
	end
	
end

