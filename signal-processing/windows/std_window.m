% 2016-03-03 21:37:35.253821897 +0100
% Karl Kastner, Berlin
%% moving block standard deviation
function s = std_moving(x,n)
	x0 = x;
	[x sdx] = sort(x);
	n = round(n/2);
	s=zeros(size(x));
	m = length(x);
	for idx=1:m
		s(idx) = std(x(max(1,idx-n):min(m,idx+n)));
	end
	x(sdx) = x;
	s(sdx) = s;
