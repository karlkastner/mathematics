% 2015-06-05 13:47:52.025537476 +0200
function X = maxfilt1(X,n)
	o = ones(n,1);
	for idx=1:size(X,2)
		X(:,idx) = ordfilt2(X(:,idx),n,o);
	end
	
end

