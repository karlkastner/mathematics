% Thu  7 Dec 13:02:26 CET 2023
function [as,flag,fls,kls] = line_search_backtracking(x0,f0,g0,fun,pk,as)
	rhols = 0.7;
	% uneducated guess
	c = 0.0001;
	% note, > sign inverted
	kls = 1;
	% TODO no magic numbers
	maxitls = 100;
%	if (pk'*g0>0)
%		error('not a descend direction');
%	end
	flag = 0;
	while (1)
		fls = fun(x0+as*pk);
		if (fls < f0 + c*as*pk'*g0)
			% convergence
			break;
		end
		as  = rhols*as;
		kls = kls+1;
		if (as < 10*eps)
			flag = -1;
			break;
		end % as < sqrt(eps)
	end % while 1
end % line_search_backtracking


