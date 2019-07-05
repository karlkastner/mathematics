% 2019-03-14 15:47:46.781449521 +0100 swapoddeven.m
% Karl Kastner, Berlin
%% swap odd and even elements in a vector
% TODO vectorize
function x = swapevenodd(x)
	for idx=1:2:size(x,1)-1
		x_ = x(idx,:);
		x(idx,:) = x(idx+1,:);
		x(idx+1,:) = x_;
	end
end

