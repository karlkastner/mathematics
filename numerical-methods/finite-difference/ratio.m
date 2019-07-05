% Mon 21 Nov 00:20:50 CET 2016
%% ratio of two subsequent values
function x = ratio(x)
	if (isvector(x))
		x = x(1:end-1)./x(2:end);
	else
		x = x(1:end-1,:)./x(2:end,:);
	end
end
