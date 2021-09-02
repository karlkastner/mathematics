% 2021-06-14 19:43:28.849599603 +0200
% TODO allow for circular
function AA=shuffle_index(A,s)
	n = size(A);
	AA = zeros(n(1),n(2));
	 for idx=1:size(A,2)
	  for jdx=1:size(A,1)
		 AA(idx,jdx) = interp2((1:n(2))',(1:n(1)),A,min(n(2),max(1,idx+s*randn())),max(1,min(n(1),jdx+s*randn())));
	  end
	 end
	A=AA.';
end

