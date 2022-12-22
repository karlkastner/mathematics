% 2018-10-21 12:21:53.205319085 +0200
%% filter columns of x (matlab does only support vector input)
function x = filter2(a,b,x)
	for idx=1:size(x,2)
		x(:,idx) = filter(a,b,x(:,idx));
	end
end

