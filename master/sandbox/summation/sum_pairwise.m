% Fri Sep 23 00:27:12 MSD 2011
% Karl Kästner, Berlin

function x = sum_pairwise(x)
	l = int32(length(x));
	while(1)
	% todo no interleaving, first half + second halfe (faster)

		l_sum = idivide(l, int32(2));
		l = l_sum + bitand(double(l),1);
		x(1:l_sum) = x(1:l_sum) + x(l+1:l+l_sum);
		if (1 == l)
			x = x(1);
			return;
		end
%{
		if (bitand(length(x),1))
			% take care of the odd element
			odd = x(end);
			x = x(1:2:end-1) + x(2:2:end);
			x(end) = x(end) + odd;
		else
			% no odd element
			x = x(1:2:end-1) + x(2:2:end);
		end
		if (1 == length(x))
			break;
		end
%}
	end
end % function sum_pairwise

