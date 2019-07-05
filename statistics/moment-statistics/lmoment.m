% Thu Dec 11 14:21:16 CET 2014
% Karl Kastner, Berlin
%
%% l-moment of vector x
function l = lmoment(x,mdx)
	% compute order statistics
	x = sort(x);
	n = length(x);
	switch(mdx)
	case {1}
		l = mean(l);
	case {2}
		id = (2:n-1)';
		%w = 0.5/nchoosek(n,2)*(nchoosek_man(idx-1,ones(n,1)) - nchoosek_man(n-idx,ones(n,1)));
		% simplification
		w  = 1/(n*(n-1))*(2*id-n-1);
		l = w'*x(id);
	end
end % lmoment


