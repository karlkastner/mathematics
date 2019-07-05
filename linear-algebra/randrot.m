% 2016-01-15 09:56:50.642149576 +0100
%% random rotation matrix
function R = randrot(n)
	R = eye(n);
	for idx=1:n
	 for jdx=idx+1:n
		alpha = 2*pi*rand();
		R = R*rot2(alpha,idx,jdx,n);
	 end
	end
end

