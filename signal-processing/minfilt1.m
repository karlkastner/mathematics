% 2015-06-05 13:47:37.345987364 +0200
function X = minfilt1(X,n)
	o = ones(n,1);
	for idx=1:size(X,2)
		X(:,idx) = ordfilt2(X(:,idx),1,o);
	end
	
end

